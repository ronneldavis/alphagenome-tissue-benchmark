import json,sys
SP=sys.argv[1]
D=json.load(open(f"{SP}/chart_data.json"))
loci=D['t2d_loci']

TIS=['ISLET','PANCREAS','LIVER','ADIPOSE','MUSCLE']
NICE={'ISLET':'Islet / beta cell','PANCREAS':'Pancreas (bulk)','LIVER':'Liver','ADIPOSE':'Fat','MUSCLE':'Skeletal muscle'}
LDL={'ISLET':(0.345,-0.00,0.80),'PANCREAS':(0.474,0.11,0.93),'LIVER':(1.130,0.41,1.99),
     'ADIPOSE':(0.167,-0.08,0.46),'MUSCLE':(0.033,-0.13,0.26)}
T2D={'ISLET':(0.178,-0.06,0.49),'PANCREAS':(0.082,-0.07,0.29),'LIVER':(0.149,-0.09,0.51),
     'ADIPOSE':(0.298,-0.08,0.85),'MUSCLE':(0.246,-0.06,0.72)}

LO,HI=-0.35,2.10
X0,X1=104.0,424.0
def sx(v): return X0+(v-LO)/(HI-LO)*(X1-X0)

def panel(data,cls):
    rows=[];y=30.0;BH,GAP=22.0,12.0
    zero=sx(0)
    for t in TIS:
        m,lo,hi=data[t]; cy=y+BH/2
        bx=min(sx(m),zero); bw=abs(sx(m)-zero)
        rows.append(f'<g class="row" tabindex="0" data-t="{NICE[t]}" data-m="{m:.3f}" data-ci="{lo:.2f} to {hi:.2f}">'
          f'<rect class="hit" x="{X0-96:.1f}" y="{y:.1f}" width="{X1-X0+96:.1f}" height="{BH:.1f}"></rect>'
          f'<text class="tlab" x="{X0-12:.1f}" y="{cy+4:.1f}">{NICE[t]}</text>'
          f'<rect class="bar {cls}" x="{bx:.1f}" y="{y+4:.1f}" width="{max(bw,1.2):.1f}" height="{BH-8:.1f}" rx="3"></rect>'
          f'<line class="ci" x1="{sx(lo):.1f}" y1="{cy:.1f}" x2="{sx(hi):.1f}" y2="{cy:.1f}"></line>'
          f'<line class="cap" x1="{sx(lo):.1f}" y1="{cy-4:.1f}" x2="{sx(lo):.1f}" y2="{cy+4:.1f}"></line>'
          f'<line class="cap" x1="{sx(hi):.1f}" y1="{cy-4:.1f}" x2="{sx(hi):.1f}" y2="{cy+4:.1f}"></line>'
          f'<text class="vlab" x="{sx(hi)+7:.1f}" y="{cy+4:.1f}">{m:.2f}</text></g>')
        y+=BH+GAP
    ticks=''.join(f'<line class="grid" x1="{sx(v):.1f}" y1="24" x2="{sx(v):.1f}" y2="{y-GAP:.1f}"></line>'
                  f'<text class="tick" x="{sx(v):.1f}" y="{y-GAP+16:.1f}">{v:g}</text>' for v in (0,0.5,1.0,1.5,2.0))
    return (f'<svg viewBox="0 0 470 {y-GAP+24:.0f}" role="img">{ticks}'
            f'<line class="zero" x1="{sx(0):.1f}" y1="24" x2="{sx(0):.1f}" y2="{y-GAP:.1f}"></line>'
            +''.join(rows)+'</svg>')

# loci chart
W,H=660,236; PL,PR,PT,PB=34,14,16,26
vals=[l['v'] for l in loci]; n=len(vals); TOP=14.0
bw=(W-PL-PR)/n
def ly(v): return PT+(1-max(v,0)/TOP)*(H-PT-PB)
bars=[]
for i,l in enumerate(loci):
    x=PL+i*bw; h=max(ly(0)-ly(l['v']),0.8)
    hot = l['v']>1
    bars.append(f'<rect class="lb{" hot" if hot else ""}" x="{x:.2f}" y="{ly(l["v"]):.2f}" width="{max(bw-1.1,1.1):.2f}" '
                f'height="{h:.2f}" tabindex="0" data-g="{l["gene"] or l["rsid"]}" data-r="{l["rsid"]}" data-v="{l["v"]:.2f}"></rect>')
tcf_i=next(i for i,l in enumerate(loci) if l['rsid']=='rs7903146')
tcf_x=PL+tcf_i*bw+bw/2
gy=[(0,'0'),(1,'1'),(5,'5'),(10,'10')]
grid=''.join(f'<line class="grid" x1="{PL}" y1="{ly(v):.1f}" x2="{W-PR}" y2="{ly(v):.1f}"></line>'
             f'<text class="tick ty" x="{PL-8}" y="{ly(v)+4:.1f}">{lab}</text>' for v,lab in gy)
loci_svg=(f'<svg viewBox="0 0 {W} {H}" role="img">{grid}'
  f'<line class="thresh" x1="{PL}" y1="{ly(1):.1f}" x2="{W-PR}" y2="{ly(1):.1f}"></line>'
  +''.join(bars)+
  f'<g class="ann"><line x1="{PL+bw*0.5:.1f}" y1="{ly(13.32)-6:.1f}" x2="{PL+bw*7:.1f}" y2="{ly(13.32)-6:.1f}"></line>'
  f'<text x="{PL+bw*7.8:.1f}" y="{ly(13.32)-2:.1f}">PLEKHA1 — but it lights up in <tspan class="em">every</tspan> tissue equally</text></g>'
  f'<g class="ann"><line x1="{PL+bw*1.5:.1f}" y1="{ly(5.21):.1f}" x2="{PL+bw*8:.1f}" y2="{ly(5.21):.1f}"></line>'
  f'<text x="{PL+bw*8.8:.1f}" y="{ly(5.21)+4:.1f}">MTNR1B — islet-specific, and a real beta-cell gene</text></g>'
  f'<g class="ann tcf"><line x1="{tcf_x:.1f}" y1="{ly(0.25):.1f}" x2="{tcf_x:.1f}" y2="{ly(0.25)-26:.1f}"></line>'
  f'<text x="{tcf_x+6:.1f}" y="{ly(0.25)-30:.1f}">TCF7L2 rs7903146 sits here</text></g>'
  f'<text class="tick" x="{PL}" y="{H-6}">118 fine-mapped diabetes signals, ranked by their strongest effect in any tissue</text>'
  '</svg>')

tab=''.join(f'<tr><th scope="row">{NICE[t]}</th><td class="n">{T2D[t][0]:.2f}</td>'
            f'<td class="n mu">{T2D[t][1]:.2f} to {T2D[t][2]:.2f}</td>'
            f'<td class="n">{LDL[t][0]:.2f}</td><td class="n mu">{LDL[t][1]:.2f} to {LDL[t][2]:.2f}</td></tr>' for t in TIS)

html=open(f"{SP}/tpl.html").read()
for k,v in [('__P_LDL__',panel(LDL,'ldl')),('__P_T2D__',panel(T2D,'t2d')),
            ('__LOCI__',loci_svg),('__TABLE__',tab),
            ('__N_ABOVE__',str(D['t2d_above1'])),('__N_T2D__',str(D['n_t2d']))]:
    html=html.replace(k,v)
open(f"{SP}/report.html","w").write(html)
print("wrote report.html", len(html), "bytes")
