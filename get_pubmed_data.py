from Bio import Entrez

pubmed_search_results = {
    "search": "",
    "result-dictionaries": []
}


def pubmed_search(search_term_pubmed='intellectual disability',
                  gen_symbols_pubmed='PRRT2'):
    combined_search = "" + search_term_pubmed
    split_gen_sym = gen_symbols_pubmed.split(",")
    for word in split_gen_sym:
        combined_search += (" AND " + word)

    pubmed_search_results["search"] = (
    "https://pubmed.ncbi.nlm.nih.gov/?term={}".format(
        combined_search.replace(" ", "+")))

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
            

            all_key_words = ', '.join(
                [key_word for key_words in article_key_words for key_word in
                 key_words])

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
                "Keywords": all_key_words,
                "Abstract": article_abstract,
                "PMID": article_pmid,
                "Link": article_link
            }

            pubmed_search_results["result-dictionaries"].append(results_dict)

        except:
            pass



    print(pubmed_search_results["search"])
    for item in range(0,len(pubmed_search_results["result-dictionaries"])):
        print(str(pubmed_search_results["result-dictionaries"][item]))



pubmed_search()
