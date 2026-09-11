import json,urllib.request,urllib.error,time,sys
URL="https://api.platform.opentargets.org/api/v4/graphql"
STUDY="GCST006867"
def q(query,v=None,tries=4):
    for t in range(tries):
        req=urllib.request.Request(URL,data=json.dumps({"query":query,"variables":v or {}}).encode(),
            headers={"Content-Type":"application/json"})
        try: return json.load(urllib.request.urlopen(req,timeout=180))
        except urllib.error.HTTPError as e:
            body=e.read().decode()[:400]
            if t==tries-1: print("ERR:",body); raise
            time.sleep(3)
        except Exception:
            if t==tries-1: raise
            time.sleep(3)

Q="""query($s:[String!], $i:Int!, $n:Int!){
 credibleSets(studyIds:$s, page:{index:$i,size:$n}){
  count
  rows{
   studyLocusId chromosome position pValueMantissa pValueExponent finemappingMethod confidence
   variant{ id chromosome position referenceAllele alternateAllele rsIds mostSevereConsequence{ id label } }
   locus(page:{index:0,size:200}){ count rows{ posteriorProbability is95CredibleSet
       variant{ id chromosome position referenceAllele alternateAllele rsIds mostSevereConsequence{ id label } } } }
   l2GPredictions(page:{index:0,size:3}){ rows{ score target{ approvedSymbol } } }
  } } }"""

allrows=[]; i=0
while True:
    r=q(Q,{"s":[STUDY],"i":i,"n":25})
    d=r["data"]["credibleSets"]
    allrows+=d["rows"]
    print(f"page {i}: got {len(d['rows'])} / total {d['count']}", flush=True)
    i+=1
    if len(allrows)>=d["count"] or not d["rows"]: break
    time.sleep(0.5)
json.dump(allrows, open(f"{sys.argv[1]}/t2d_credsets.json","w"))
print("saved", len(allrows), "credible sets")
