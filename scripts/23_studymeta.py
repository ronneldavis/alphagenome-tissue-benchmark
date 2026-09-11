import json,urllib.request,sys,time
URL="https://api.platform.opentargets.org/api/v4/graphql"
def q(s):
    r=urllib.request.Request(URL,data=json.dumps({"query":s}).encode(),headers={"Content-Type":"application/json"})
    return json.load(urllib.request.urlopen(r,timeout=120))
panel=json.load(open(f"{sys.argv[1]}/trait_panel.json"))
ids=[p['study'] for p in panel]+["GCST006867"]
out={}
for sid in ids:
    try:
        d=q('{ study(studyId:"%s"){ id publicationFirstAuthor publicationDate publicationJournal pubmedId nSamples nCases traitFromSource projectId } }'%sid)["data"]["study"]
        out[sid]=d
        print(f"{sid:<20s} {str(d.get('publicationFirstAuthor')):<16s} {str(d.get('publicationDate'))[:4]:<6s} pmid={d.get('pubmedId')} n={d.get('nSamples')} | {(d.get('traitFromSource') or '')[:38]}")
    except Exception as e: print(sid,"ERR",str(e)[:60])
    time.sleep(0.2)
json.dump(out,open(f"{sys.argv[1]}/study_meta.json","w"))
