import itertools
import re
import os

from nltk.tokenize import sent_tokenize


class CoOccurrence:
    """
    A class to calculate the co-occurrence of titles and abstracts of articles.
    It creates combinations of gene symbols and phentopyes which will be
    used to search for new relationsips in articles. It takes about 7 minutes
    to calculate all of the possible (41709140) combinations.
    :param data: a list that contains dictionary with the info of an article.
    :param gene_symbols: the file that contains the gene symbols.
    :param phenotypes: the file that contains the phenotypes.
    :param title: boolean whether or not to include the title for the
    calculation. Default is true.
    :param in_sentence: boolean whether or not to inculde a sentence of the
    abstract for the calculation. Default is true.
    :param in_abstract: boolean whether or not to inculde the abstract for
    the calculation. Default is true.
    :param multiple_times_in_abstract: boolean whether or not to inculde if a
    combination occurs in multiple abstracts for the calculation. Default is
    true.
    """

    def __init__(self, data, title=True,
                 in_sentence=True, in_abstract=True,
                 multiple_times_in_abstract=True):
        self.__directory = f'co_occurrence_algorithm{os.sep}data_files'
        self.data_list = data
        self.gene_symbols = os.path.join(self.__directory,'GenPanels_merged_DG-2.17.0.tsv')
        self.phenotypes = os.path.join(self.__directory, 'gene_condition_source_id.txt')
        self.combinations = []
        self.title = title
        self.in_sentence = in_sentence
        self.in_abstract = in_abstract
        self.multiple_times_in_abstract = multiple_times_in_abstract
        self.tokenized = []
        self.co_occurrence = {}
        self.__get_combinations()

    def pre_process_data(self):
        """
        Preprocessing of the data. The abstract is split up into sentences via
        the sent_tokenize funtion of the nltk toolkit. It creates a new
        dictionary with the preprocesses titles and abstracts.
        """
        for data in self.data_list:
            data_dict = {}
            title = data['title']
            data_dict['title'] = title.strip()
            abstract = self.__pre_process_abstract(data['abstract']).strip()
            abstract_sentences = sent_tokenize(abstract)
            data_dict['abstract'] = abstract_sentences
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
            for combination in self.combinations:
                gene_space = f' {combination[0]} '
                gene_parentheses = f'({combination[0]})'
                phenotype = combination[1].lower()
                if self.title:
                    self.__check_part(gene_space, gene_parentheses,
                                      phenotype, combination, data['title'],
                                      10)

                if self.in_sentence:
                    for sentence in data['abstract']:
                        self.__check_part(gene_space, gene_parentheses,
                                          phenotype, combination, sentence, 5)
                if self.in_abstract:
                    abstract = ' '.join(data['abstract'])
                    self.__check_part(gene_space, gene_parentheses,
                                      phenotype, combination, abstract, 2)

        if self.multiple_times_in_abstract:
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
                    try:
                        self.co_occurrence[combination] += 1
                    except KeyError:
                        self.co_occurrence[combination] = 1

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
                     combination, part, amount):
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
        :param combination: A combination of a gene symbol and phenotype,
        e.g.: (PRRT2, Intelectual disability)
        :param part: A part of the article, which can be either the title,
        abstract or a sentence of an abstract.
        :param amount: An amount of points which will be given when there is a
        match. Title: 10, sentence: 5, abstract: 2, multiple abstracts: 1.
        """
        if (gene_space in part
            or gene_parentheses in part) \
                and phenotype in part:
            try:
                self.co_occurrence[combination] += amount
            except KeyError:
                self.co_occurrence[combination] = amount

    def __get_combinations(self):
        """
        Reads a file that contains the HGNC gene symbols and adds that to a
        list. Reads a file that contains phenotypes and adds that to a list.
        Then all of the possible combinations will be generated and saved to
        a new list.
        """
        genes = []
        phenotypes = []
        with open(self.gene_symbols, 'r') as genes_file:
            genes_file.readline()
            for line in genes_file:
                gene = line.split('\t')[0]
                genes.append(gene.strip())
        with open(self.phenotypes, 'r') as phenotype_file:
            phenotype_file.readline()
            for line in phenotype_file:
                phenotype = line.split('\t')[4]
                phenotypes.append(phenotype.strip())

        self.combinations = list(itertools.product(genes, phenotypes))
