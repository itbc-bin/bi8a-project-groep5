from Bio import Entrez
from flask import Flask, render_template, request

app = Flask(__name__)


# TODO
# return dictionary met search_term+gensymbool
# return lijst met alle sites ( searchterm+gensymbool)
# return
# functie om verbinding te maken met de database ( local database )

@app.route('/')
def hello_world():
    return render_template('index.html')


@app.route('/test')
def pubmed_search():
    search_term_pubmed = request.args.get('search_term')
    if search_term_pubmed is None: search_term_pubmed = "shark fish green"  # axolotl
    html_code = """
                <a href=/> HOME</a><br>
                <form method="get">
                <input type="text" name="search_term">
                <textarea name="gen_symbols"></textarea>
                <input type="submit" value="Submit">
                </form>
                <hr>
                """

    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='10000',
                            retmode='xml',
                            term=search_term_pubmed)
    results_handle = Entrez.read(handle)
    ids = ','.join(results_handle['IdList'])
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results_fetch = Entrez.read(handle)
    paper = results_fetch

    # {'titel': 'de titel', 'auteur': 'de auteur', 'publicatie_datum': 'datum',
    # 'keywords': 'de keywords', 'abstract': 'de abstract', 'pmid': '31801583', 'link': 'de link'}

    for i, paper in enumerate(paper['PubmedArticle']):

        article_title = paper['MedlineCitation']['Article']['ArticleTitle']
        article_authors = []
        for i in paper['MedlineCitation']['Article']['AuthorList']:
            last_name = i['LastName']
            fore_name = i['ForeName']
            initials = i['Initials']
            author = "{} {} {}".format(last_name,
                                       fore_name,
                                       initials)
            article_authors.append(author)

        article_pub_year = \
            paper['MedlineCitation']['Article']['Journal']['JournalIssue'][
                'PubDate']['Year']
        article_key_words = paper['MedlineCitation']['KeywordList']
        try:
            article_abstract = paper['MedlineCitation']['Article']['Abstract']
        except:
            article_abstract = "-"
        article_pmid = paper['MedlineCitation']['PMID']
        article_link = "-"

        results_dict = {
            "Title": article_title,
            "Authors": article_authors,
            "Publication_year": article_pub_year,
            "Keywords": article_key_words,
            "Abstract": article_abstract,
            "PMID": article_pmid,
            "Link": article_link
        }
        print(results_dict)
        print("\n" * 2)

    return html_code


if __name__ == '__main__':
    app.run()
