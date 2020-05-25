import warnings
import os

import mysql.connector
from Bio import Entrez

warnings.filterwarnings('ignore', message='Numerical issues were encountered')


class PubmedSearch:

    def __init__(self, e_mail, search_word, gene_symbols, year=2000):
        self.email = e_mail
        self.search_word = search_word
        self.gene_symbols = gene_symbols
        self.year = year
        self.articles_data = []
        self.articles = []
        self.results = []

    def search_pubmed(self):
        Entrez.email = self.email
        for gene_symbol in self.gene_symbols.split(os.linesep):
            gene_symbols = self.get_symbols_from_database(gene_symbol)
            all_aliases = gene_symbols.split(',')
            search_term = f'{self.search_word} AND ({gene_symbol}'
            for counter, symbol in enumerate(all_aliases):
                if counter == 0:
                    continue
                if counter == len(all_aliases) - 1:
                    search_term += f' OR {symbol.strip()})'
                    break
                search_term += f' OR {symbol.strip()}'
            search_term += f' AND ("{self.year}"[Date - Publication] : ' \
                           f'"3000"[Date - Publication])'
            term_handle = Entrez.esearch(db='pubmed',
                                         sort='relevance',
                                         retmax='10000',
                                         retmode='xml',
                                         term=search_term)

            ids = Entrez.read(term_handle)['IdList']
            amount_ids = len(ids)
            search_link = search_term.replace(' ', '+') \
                .replace('(', '%28').replace(')', '%29')
            link_page = f'https://pubmed.ncbi.nlm.nih.gov/?term={search_link}'
            self.articles_data.append(
                {'zoekwoord': gene_symbol, 'gezochte_symbolen': gene_symbols,
                 'aantal_hits': amount_ids, 'link': link_page})
            if ids:
                id_handle = Entrez.efetch(db='pubmed',
                                          retmode='xml',
                                          id=ids)
                articles = Entrez.read(id_handle)['PubmedArticle']
                self.articles.append({gene_symbol: articles})

    def parse_results(self):
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
            cursor.execute(
                "insert into articles(pubmed_id, title, publication_year, "
                "keywords, article_link, abstract, authors) values ('{}',"
                "'{}','{}','{}','{}','{}','{}')".format(pmid,
                                                        title, pub_year,
                                                        keywords, link,
                                                        abstract,
                                                        authors))
            connection.commit()
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
        return ', '.join(results[0])


if __name__ == '__main__':
    term = 'intellectual disability'
    words = 'PRRT2\nKCNMA1'
    email = 'christiaanposthuma@gmail.com'
    search = PubmedSearch(e_mail=email, search_word=term, gene_symbols=words,
                          year=2015)
    search.search_pubmed()
    search.parse_results()
    search.insert_to_database()