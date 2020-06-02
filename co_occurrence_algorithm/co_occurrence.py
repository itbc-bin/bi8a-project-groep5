import itertools
import os
import re
import time

import mysql.connector
import nltk
from nltk.tokenize import sent_tokenize

NLTK_DIR = os.path.join(os.getcwd(), 'nltk_data')
nltk.data.path.append(NLTK_DIR)

"""
Class to perform the co occurrence algorithm.
@author: Yaris
"""


class CoOccurrence:
    """
    Class to perform the co occurrence algorithm. The algorithm will search
    combinations in all of the articles which were found with the PubMed
    search. The algorithm generates combninations of the phenotype (from user
    input) and 4409 different gene symbols from the GenPanels file. The
    combinations are scored on different aspects. If the combination occurs in
    a title, that combination will get 10 points, if the combination occurs in
    a single sentence of the abstract, the combination will get 5 points. If
    the combination occurs in the abstract, the combination will get 2 points.
    Lastly, if the combination occurs in multiple abstracts, the combination
    will get 1 point.
    """

    def __init__(self, data, url_id, term, title=None, in_title=True,
                 in_sentence=True, in_abstract=True,
                 in_multiple_abstracts=True):
        """
        Initialization method of the co occurrence class.
        :param data: The data containing the articles.
        :param url_id: The id of the url.
        :param term: The search term from user input.
        :param title: The title for the search from user input.
        :param in_title: Boolean whether or not to include the title for the
        calculation. Default is true.
        :param in_sentence: Boolean whether or not to inculde a sentence of the
        abstract for the calculation. Default is true.
        :param in_abstract: Boolean whether or not to inculde the abstract for
        the calculation. Default is true.
        :param in_multiple_abstracts: Boolean whether or not to inculde if a
        combination occurs in multiple abstracts for the calculation. Default is
        true.
        """
        self.__directory = os.path.join('co_occurrence_algorithm',
                                        'data_files')
        self.data_list = data
        self.url_id = url_id
        self.phenotype = [term]
        self.gene_symbols = os.path.join(self.__directory,
                                         'GenPanels_merged_DG-2.17.0.tsv')
        self.combinations = []
        self.title = title
        self.in_title = in_title
        self.in_sentence = in_sentence
        self.in_abstract = in_abstract
        self.in_multiple_abstracts = in_multiple_abstracts
        self.options = [{'in title': in_title, 'in sentence': in_sentence,
                         'in abstract': in_abstract,
                         'multiple times in abstract': in_multiple_abstracts}]
        self.tokenized = []
        self.co_occurrence = {}
        self.__download_punkt()
        self.__get_combinations()

    def pre_process_data(self):
        """
        Preprocessing of the data. The abstract is split up into sentences via
        the sent_tokenize funtion of the nltk toolkit. It creates a new
        dictionary with the preprocesses titles and abstracts.
        """
        for data in self.data_list:
            data_dict = {}
            title = data['Title']
            data_dict['title'] = title.strip()
            abstract = self.__pre_process_abstract(data['Abstract']).strip()
            abstract_sentences = sent_tokenize(abstract)
            data_dict['abstract'] = abstract_sentences
            data_dict['pmid'] = data['PMID']
            self.tokenized.append(data_dict)

    def calculate_co_occurrence(self):
        """
        The base algorithm for co-occurrence. It calculates 4 different scores,
        based on co-occurence. A combination consists of a gene symbol, which
        has two options: 1.  PRRT2  or 2. (PRRT2), and a phenotype: 3
        e.g. intellectual disability.
        The user has the option to choose what will be scored. There is an
        option to give 10 points if a combination occurs in the title, an
        option to give 5 points if a combination occurs in a sentence, an
        option to give 2 points if a combination occurs in an abstract and
        finally an option to give 1 point if a combination occurs in multiple
        different abstracts. A dictionary will be made that keeps track of
        combinations with a score higher than 0.
        """
        for data in self.tokenized:
            abstract = ' '.join(data['abstract'])
            for combination in self.combinations:
                gene_space = f' {combination[0]} '
                gene_parentheses = f'({combination[0]})'
                phenotype = combination[1].lower()
                if self.in_title:
                    self.__check_part(gene_space, gene_parentheses,
                                      phenotype, data['title'],
                                      10, data['pmid'])

                if self.in_sentence:
                    for sentence in data['abstract']:
                        self.__check_part(gene_space, gene_parentheses,
                                          phenotype, sentence, 5, data['pmid'])
                if self.in_abstract:
                    self.__check_part(gene_space, gene_parentheses,
                                      phenotype, abstract, 2, data['pmid'])

        if self.in_multiple_abstracts:
            all_abstracts = [' '.join(data['abstract']) for data in
                             self.tokenized]
            for combination in self.combinations:
                amount = 0
                gene_space = f' {combination[0]} '
                gene_parentheses = f'({combination[0]})'
                phenotype = combination[1].lower()
                for abstract in all_abstracts:
                    if (gene_space in abstract
                        or gene_parentheses in abstract) \
                            and phenotype in abstract:
                        amount += 1
                if amount > 1:
                    gene = gene_space.strip()
                    if f'{gene}, {phenotype}' not in self.co_occurrence:
                        self.co_occurrence[f'{gene}, {phenotype}'] = {}
                        self.co_occurrence[f'{gene}, {phenotype}'][
                            'options'] = self.options
                    try:
                        self.co_occurrence[f'{gene}, {phenotype}'][
                            'amount'] += 1
                    except KeyError:
                        self.co_occurrence[f'{gene}, {phenotype}'][
                            'amount'] = 1

    def save_to_db(self):
        """
        Save co occurrence results to a database, so the results can always
        be retrieved.
        """
        connection = self.__connection_database()
        cursor = connection.cursor()
        cursor.execute("select count(*) from algorithm_results")
        result_id = cursor.fetchone()[0] + 1
        cursor.execute(
            "insert into results (url_id, title) values ('{}', '{}')".format(
                self.url_id, self.title))
        connection.commit()

        cursor.execute(
            "select id from results where url_id = '{}'".format(self.url_id))
        url_id_db = cursor.fetchone()[0]
        for key, value in self.co_occurrence.items():
            cursor.execute("insert into algorithm_results "
                           "(results_id, combination, amount) values ('{}', "
                           "'{}', '{}')".format(url_id_db, key,
                                                value['amount']))
            connection.commit()
            for pmid in value['ids']:
                cursor.execute(
                    "insert into pmid_result (algorithm_results_id, "
                    "pmid) values ('{}', '{}')".format(result_id, pmid))
                connection.commit()
            result_id += 1
        connection.close()

    def get_co_occurence(self):
        """
        Get the calculated co-occurence dictionary.
        :return: Calculated co-occurence dictionary.
        """
        return self.co_occurrence

    @staticmethod
    def __pre_process_abstract(abstract_data):
        """
        Pre processing of the abstract. A regex pattern is used to search for
        gene symbols, so that the gene symbols will be uppercase and the rest
        of the abstract will be lowercase.
        :param abstract_data:
        :return: a pre processed abstract.
        """
        gene_symbol_pattern = r'((?:|\()(?:[A-Z](?:[A-Z]|[0-9]|-)+)(?:|\)|,))'
        gene_symbols = re.findall(gene_symbol_pattern, abstract_data)
        abstract_data = abstract_data.lower()
        for gene_symbol in gene_symbols:
            abstract_data = abstract_data.replace(gene_symbol.lower(),
                                                  gene_symbol)
        return abstract_data

    def __check_part(self, gene_space, gene_parentheses, phenotype,
                     part, amount, pmid):
        """
        Checks is either one of the gene symbol types and the phenotype is
        present in a part of the article, which can be the title, abstract or
        a sentence of the abstract. The co-occurrence dictionary will be
        updated with the corresponding value when there is a match.
        :param gene_space: A gene symbol which has spaces around it,
        e.g.:  PRRT2 .
        :param gene_parentheses: A gene symbol that has parentheses around it,
        e.g.: (PRRT2).
        :param phenotype: A phenotype, e.g.: intellectual disability.
        :param part: A part of the article, which can be either the title,
        abstract or a sentence of an abstract.
        :param amount: An amount of points which will be given when there is a
        match. Title: 10, sentence: 5, abstract: 2, multiple abstracts: 1.
        """
        if (gene_space in part
            or gene_parentheses in part) \
                and phenotype in part:
            gene = gene_space.strip()
            if f'{gene}, {phenotype}' not in self.co_occurrence:
                self.co_occurrence[f'{gene}, {phenotype}'] = {}
                self.co_occurrence[f'{gene}, {phenotype}'][
                    'options'] = self.options
            try:
                self.co_occurrence[f'{gene}, {phenotype}']['amount'] += amount
                if pmid not in self.co_occurrence[f'{gene}, {phenotype}'][
                    'ids']:
                    self.co_occurrence[f'{gene}, {phenotype}']['ids'].append(
                        pmid)
            except KeyError:
                self.co_occurrence[f'{gene}, {phenotype}']['amount'] = amount
                self.co_occurrence[f'{gene}, {phenotype}']['ids'] = [pmid]

    def __get_combinations(self):
        """
        Reads a file that contains the HGNC gene symbols and adds that to a
        list. Reads a file that contains phenotypes and adds that to a list.
        Then all of the possible combinations will be generated and saved to
        a new list.
        """
        connection = self.__connection_database()
        cursor = connection.cursor()
        cursor.execute("select gene_symbol from gene_symbols")
        gene_symbols = cursor.fetchall()
        genes = [gene_symbol[0].strip() for gene_symbol in gene_symbols]
        self.combinations = list(itertools.product(genes, self.phenotype))

    @staticmethod
    def __connection_database():
        """
        Make connection to the database
        :return: The connection
        """
        connection = mysql.connector.connect(
            host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
            user='owe7_pg5@hannl-hlo-bioinformatica-mysqlsrv',
            database='owe7_pg5',
            password='blaat1234')

        return connection

    @staticmethod
    def __download_punkt():
        """
        Download the punkt directory from nltk if it has not been installed
        yet.
        """
        nltk_dir = os.path.join(os.getcwd(), 'nltk_data')

        if not os.path.isdir(nltk_dir):
            nltk.download('punkt', download_dir=nltk_dir)
            time.sleep(15)
