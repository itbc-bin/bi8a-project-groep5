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
    gen_symbols_pubmed = request.args.get('gen_symbols')
    if search_term_pubmed is None: search_term_pubmed = "intellectual disability"
    if gen_symbols_pubmed is None: gen_symbols_pubmed = "PRRT2"

    combined_search = "" + search_term_pubmed
    split_gen_sym = gen_symbols_pubmed.split(",")
    for word in split_gen_sym:
        combined_search += (" AND " + word)


    html_code = """
                <a href=/> HOME</a><br>
                <form method="get">
                Search term: <input type="text" name="search_term"><br>
                Gene symbols: <input type="text" name="gen_symbols">
                <input type="submit" value="Submit">
                </form>
                
                For gen_symbols: split the terms by comma
                <hr>
                """

    Entrez.email = 'christiaanposthuma@gmail.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='10000',
                            retmode='xml',
                            term=combined_search)
    results_handle = Entrez.read(handle)
    ids = ','.join(results_handle['IdList'])
    Entrez.email = 'christiaanposthuma@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results_fetch = Entrez.read(handle)
    paper = results_fetch

    for i, paper in enumerate(paper['PubmedArticle']):
        try:
            article_title = paper['MedlineCitation']['Article']['ArticleTitle']
            article_authors = []

            for i in paper['MedlineCitation']['Article']['AuthorList']:
                try:
                    last_name = i['LastName']
                    fore_name = i['ForeName']
                    initials = i['Initials']
                    author = "{} {} {}".format(last_name,
                                               fore_name,
                                               initials)
                except:
                    author = "-"

                article_authors.append(author)

            article_pub_year = \
                paper['MedlineCitation']['Article']['Journal']['JournalIssue'][
                    'PubDate']['Year']
            article_key_words = paper['MedlineCitation']['KeywordList']

            try:
                abstract_info = \
                paper['MedlineCitation']['Article']['Abstract'][
                    'AbstractText']
                article_abstract = ' '.join(abstract_info)
            except:
                article_abstract = '--'

            article_pmid = paper['MedlineCitation']['PMID'].split(",")[0]
            article_link = "https://pubmed.ncbi.nlm.nih.gov/20190110/"

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
            # https://pubmed.ncbi.nlm.nih.gov/?term=disability+AND+autism+AND+baby

        except:
            pass

    overall_link = ("https://pubmed.ncbi.nlm.nih.gov/?term={}".format(combined_search.replace(" ","+")))
    print(overall_link)

    # return results_dict,overall_link
    return html_code


if __name__ == '__main__':
    app.run()
