import warnings
from datetime import datetime

import mysql.connector
from Bio import Entrez

warnings.filterwarnings('ignore', message='Numerical issues were encountered')

# TODO write documentation
# TODO handle bugs

class PubmedSearch:

    def __init__(self, e_mail, search_word, gene_symbols,
                 date='January 1, 2000'):
        self.email = e_mail
        self.api_key = '70603012ca2859e88695f0dae2d6dc988308'
        self.search_word = search_word.strip()
        self.gene_symbols = gene_symbols
        self.date = datetime.strptime(date, '%B %d, %Y').strftime('%Y/%m/%d')
        self.ids_data = []
        self.articles_data = []
        self.articles = []
        self.results = []

    def search_pubmed(self):
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        for gene_symbol in self.gene_symbols.split('\n'):
            gene_symbol = gene_symbol.strip()
            gene_symbols = self.get_symbols_from_database(gene_symbol)
            all_aliases = gene_symbols.split(',')
            if len(all_aliases) == 1:
                search_term = f'{self.search_word} AND {gene_symbol}'
            else:
                search_term = f'{self.search_word} AND ({gene_symbol}'
                for counter, symbol in enumerate(all_aliases):
                    if counter == 0:
                        continue
                    if counter == len(all_aliases) - 1:
                        search_term += f' OR {symbol.strip()})'
                        break
                    search_term += f' OR {symbol.strip()}'
            search_term += f' AND ("{self.date}"[Date - Publication] : ' \
                           f'"3000"[Date - Publication])'
            term_handle = Entrez.esearch(db='pubmed',
                                         sort='relevance',
                                         retmax='10000',
                                         retmode='xml',
                                         term=search_term)

            ids = Entrez.read(term_handle)['IdList']
            amount_ids = len(ids)
            link = search_term.replace(' ', '+').replace('(', '%28').replace(
                ')', '%29').replace('"', '%22').replace('[', '%5B').replace(
                ']', '%5D')
            url = f'https://pubmed.ncbi.nlm.nih.gov/?term={link}'
            self.articles_data.append(
                {'zoekwoord': gene_symbol.strip(),
                 'gezochte_symbolen': gene_symbols,
                 'aantal_hits': amount_ids, 'link': url})

            self.ids_data.append({'ids': ids, 'gene_symbol': gene_symbol})

    def __parse_ids(self):
        for data in self.ids_data:
            if data['ids']:
                id_handle = Entrez.efetch(db='pubmed',
                                          retmode='xml',
                                          id=data['ids'])
                articles = Entrez.read(id_handle)['PubmedArticle']
                self.articles.append({data['gene_symbol']: articles})

    def parse_results(self):
        self.__parse_ids()
        for articles in self.articles:
            for symbol, article in articles.items():
                for article_info in article:
                    try:
                        article_title = \
                            article_info['MedlineCitation']['Article'][
                                'ArticleTitle']
                        article_authors = []

                        for author in \
                                article_info['MedlineCitation']['Article'][
                                    'AuthorList']:
                            try:
                                last_name = author['LastName']
                                fore_name = author['ForeName']
                                inits = author['Initials']
                                auth_info = f'{last_name} {fore_name} {inits}'
                                article_authors.append(auth_info)

                            except KeyError:
                                article_authors.append('')
                        authors = ', '.join(article_authors)
                        article_pub_year = \
                            article_info['MedlineCitation']['Article'][
                                'Journal'][
                                'JournalIssue'][
                                'PubDate']['Year']
                        article_key_words = article_info['MedlineCitation'][
                            'KeywordList']
                        all_key_words = ', '.join(
                            [key_word for key_words in article_key_words for
                             key_word in key_words])

                        try:
                            abstract_info = \
                                article_info['MedlineCitation']['Article'][
                                    'Abstract'][
                                    'AbstractText']
                            article_abstract = ' '.join(abstract_info)
                        except KeyError:
                            article_abstract = '--'

                        article_pmid = \
                            article_info['MedlineCitation']['PMID'].split(',')[
                                0]
                        article_link = f'https://pubmed.ncbi.nlm.nih.gov/' \
                                       f'{article_pmid}/'

                        results_dict = {
                            'Title': article_title,
                            'Authors': authors,
                            'Publication_year': article_pub_year,
                            'Keywords': all_key_words,
                            'Abstract': article_abstract,
                            'PMID': article_pmid,
                            'Link': article_link
                        }
                        self.results.append(results_dict)

                    except KeyError:
                        pass

    def get_results(self):
        return self.results

    def insert_to_database(self):
        connection = self.connection_database()
        cursor = connection.cursor()
        for article in self.results:
            pmid = int(article['PMID'])
            title = article['Title'].replace("'", "''")
            pub_year = int(article['Publication_year'])
            keywords = article['Keywords'].replace("'", "''")
            link = article['Link']
            abstract = article['Abstract'].replace("'", "''")
            authors = article['Authors'].replace("'", "''")
            try:
                cursor.execute(
                    "insert into articles(pubmed_id, title, publication_year, "
                    "keywords, article_link, abstract, authors) values ('{}',"
                    "'{}','{}','{}','{}','{}','{}')".format(pmid,
                                                            title, pub_year,
                                                            keywords, link,
                                                            abstract,
                                                            authors))
                connection.commit()
            except mysql.connector.errors.IntegrityError:
                pass
        connection.close()

    @staticmethod
    def connection_database():
        """Make connection to the database
        :return: The connection
        """
        connection = mysql.connector.connect(
            host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
            user='owe7_pg5@hannl-hlo-bioinformatica-mysqlsrv',
            database='owe7_pg5',
            password='blaat1234')

        return connection

    @staticmethod
    def get_symbols_from_database(gene_symbol):
        gene_symbol = gene_symbol.strip()
        connection = PubmedSearch.connection_database()
        cursor = connection.cursor()
        cursor.execute("select gene_symbol, previous_gene_symbols, "
                       "alias_gene_symbols from gene_symbols where gene_symbol"
                       " = '{}';".format(gene_symbol))
        results = cursor.fetchall()
        connection.close()
        if results:
            return ', '.join(filter(None, results[0]))
        return gene_symbol
