import sys,re,json,numpy as np,pandas as pd
SP=sys.argv[1]
old=open(f"{SP}/tpl.html").read()
style=old[old.index("<style>"):old.index("</style>")+8]
fonts=old[old.index("<link rel"):old.index(">",old.index("<link rel"))+1]
R=pd.read_csv(f"{SP}/bench_final.csv").rename(columns={"Unnamed: 0":"trait"}).sort_values(['rank','p_exp'])
NICE={'ISLET':'Islet','LIVER':'Liver','ADIPOSE':'Fat','MUSCLE':'Muscle','BRAIN':'Brain','IMMUNE':'Immune',
 'HEART':'Heart','ARTERY':'Artery','LUNG':'Lung','INTESTINE':'Intestine','SKIN':'Skin'}

# --- chart 1: rank of expected tissue per trait (lollipop) ---
W,H=660,470; L,T,B=230,18,34
rows=list(R.itertuples()); n=len(rows); rh=(H-T-B)/n
def rx(r): return L+(r-1)/10*(W-L-58)
bars=[]
for i,r in enumerate(rows):
    cy=T+i*rh+rh/2
    sig=r.q_exp<0.05
    bars.append(
      f'<g class="row" tabindex="0" data-t="{r.trait}" data-m="rank {int(r._4 if False else r.rank)} of 11 · z={r.exp_z:.2f}" data-ci="expected: {NICE[r.expected]}">'
      f'<rect class="hit" x="0" y="{T+i*rh:.1f}" width="{W}" height="{rh:.1f}"></rect>'
      f'<text class="tlab" x="{L-14:.1f}" y="{cy+3.6:.1f}">{r.trait[:34]}</text>'
      f'<line class="ci" x1="{rx(1):.1f}" y1="{cy:.1f}" x2="{rx(r.rank):.1f}" y2="{cy:.1f}"></line>'
      f'<circle cx="{rx(r.rank):.1f}" cy="{cy:.1f}" r="4.6" fill="{"var(--teal)" if sig else "var(--null)"}"></circle>'
      f'<text class="vlab" x="{rx(r.rank)+9:.1f}" y="{cy+3.6:.1f}">{NICE[r.expected]}</text></g>')
ticks=''.join(f'<line class="grid" x1="{rx(v):.1f}" y1="{T-4}" x2="{rx(v):.1f}" y2="{H-B+2:.1f}"></line>'
              f'<text class="tick" x="{rx(v):.1f}" y="{H-B+16:.1f}">{v}</text>' for v in (1,3,5,7,9,11))
c1=(f'<svg viewBox="0 0 {W} {H}" role="img">{ticks}'+''.join(bars)+
    f'<text class="tick" x="{L}" y="{H-4}" style="text-anchor:start">position of the expected tissue among 11 &rarr;</text></svg>')

# --- chart 2: two scatters, assay depth vs recovery, visibility vs recovery ---
def scatter(xs,ys,xlab,cols,w=310,h=210):
    pl,pr,pt,pb=42,14,14,34
    x0,x1=min(xs),max(xs); y0,y1=min(ys+[0]),max(ys)
    x1=x1 if x1>x0 else x0+1
    sx=lambda v: pl+(v-x0)/(x1-x0)*(w-pl-pr)
    sy=lambda v: h-pb-(v-y0)/((y1-y0) or 1)*(h-pt-pb)
    pts=''.join(f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="4" fill="{c}" opacity=".85"></circle>'
                for a,b,c in zip(xs,ys,cols))
    gy=''.join(f'<line class="grid" x1="{pl}" y1="{sy(v):.1f}" x2="{w-pr}" y2="{sy(v):.1f}"></line>'
               f'<text class="tick ty" x="{pl-7}" y="{sy(v)+3.5:.1f}">{v:g}</text>' for v in (0,0.5,1.0,1.5))
    return (f'<svg viewBox="0 0 {w} {h}" role="img">{gy}<line class="zero" x1="{pl}" y1="{sy(0):.1f}" '
            f'x2="{w-pr}" y2="{sy(0):.1f}"></line>{pts}'
            f'<text class="tick" x="{(pl+w-pr)/2:.0f}" y="{h-6}">{xlab}</text></svg>')
cols=['var(--teal)' if q<0.05 else 'var(--null)' for q in R.q_exp]
c2=scatter(list(R.n_tracks_exp),list(R.exp_z),'tracks available for the expected tissue',cols)
c3=scatter([round(v,3) for v in R.visibility],list(R.exp_z),'does the model see the variants at all?',cols)

tab=''.join(f'<tr><th scope="row">{r.trait}</th><td>{NICE[r.expected]}</td><td class="n">{int(r.rank)}</td>'
            f'<td class="n">{r.exp_z:+.2f}</td><td class="n mu">{r.q_exp:.3f}</td></tr>' for r in rows)
open(f"{SP}/art_parts.json","w").write(json.dumps(dict(c1=c1,c2=c2,c3=c3,tab=tab,style=style,fonts=fonts)))
print("parts built; traits:",n,"| significant:",int((R.q_exp<0.05).sum()))
