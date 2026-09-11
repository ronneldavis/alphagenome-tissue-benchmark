import json,sys
SP=sys.argv[1]; P=json.load(open(f"{SP}/art_parts.json"))
body = """
<div class="wrap">
<header class="col">
  <p class="eyebrow">AlphaGenome benchmark &middot; 20 traits &middot; 11 tissues &middot; 09 Sep 2026</p>
  <h1>When is a tissue prediction trustworthy?</h1>
  <p class="lede">AI models can now read a stretch of DNA and say what it does in each tissue of the
  body. That lets you ask, of any disease variant, <strong>where in the body it acts</strong>.
  I tested how far you can trust that. Across 20 traits whose answer is already known, the model
  got the right tissue first <strong>half the time</strong> &mdash; far better than chance. The
  useful part is what separates the wins from the misses, and it isn't what I expected.</p>
  <dl class="runbar">
    <div><dt>Traits tested</dt><dd>20</dd></div>
    <div><dt>Tissues compared</dt><dd>11</dd></div>
    <div><dt>Variants scored</dt><dd>2,614</dd></div>
    <div><dt>Right tissue 1st</dt><dd>10 / 20</dd></div>
    <div><dt>Right tissue enriched</dt><dd>12 / 20</dd></div>
  </dl>
</header>

<section class="col stack">
  <h2>A correction to what I told you earlier</h2>
  <p class="callout">My first pass used one diabetes dataset and concluded the model returns a flat
  line for diabetes, probably because the pancreas is poorly measured. <strong>Both halves of that
  were wrong.</strong> With a larger diabetes study the islet ranks <em>first</em>, and the
  "poorly measured tissue" explanation fails a direct test. This page is the corrected version.</p>
</section>

<section class="col stack">
  <h2>The setup, in plain terms</h2>
  <p>Pick 20 traits where biologists already agree which tissue matters &mdash; liver enzymes act in
  the liver, heart rate in the heart, and so on. Feed the model the DNA variants linked to each
  trait. Ask it to rank all 11 tissues. Then check: did the known answer come out on top?</p>
  <p>Everything is measured against 318 random harmless mutations, so a score of
  <span class="mono">0</span> means "no different from random."</p>
</section>

<figure class="fig" style="margin-top:8px">
  <figcaption><span class="tag">Result &middot; 20 traits</span>
  <h3>Where the known answer landed</h3>
  <p class="sub">Each dot is one trait. Further left is better &mdash; position 1 means the model put
  the correct tissue first. Teal dots are traits where the correct tissue was significantly
  enriched; grey dots are traits where the model saw essentially nothing.</p></figcaption>
  __C1__
</figure>

<section class="col stack">
  <h2>Liver, heart and immune traits work. Brain traits don't.</h2>
  <p>Liver enzymes point at the liver. Heart rate points at the heart. Platelet count and Crohn's
  point at immune tissue. That's real biology recovered from DNA sequence alone, with no hints.</p>
  <p>But body mass index and cognitive ability both point at the brain in reality, and the model
  finds nothing for either. That failure is the clue.</p>
</section>

<section class="col stack">
  <h2>My explanation was wrong, and I can show it</h2>
  <p>I originally guessed that failures happen when a tissue is poorly measured &mdash; the model
  can't find what it was never shown. That's testable: if true, tissues with more data should score
  better. <strong>They don't.</strong></p>
</section>

<div class="figs">
  <figure class="fig">
    <figcaption><span class="tag"><span class="swatch" style="background:var(--null)"></span>Hypothesis 1 &mdash; rejected</span>
    <h3>Is it how well the tissue is measured?</h3>
    <p class="sub">No relationship. More data does not mean a better answer.</p></figcaption>
    __C2__
  </figure>
  <figure class="fig">
    <figcaption><span class="tag"><span class="swatch sw-ldl"></span>Hypothesis 2 &mdash; supported</span>
    <h3>Or whether it sees the variants at all?</h3>
    <p class="sub">Strong relationship. This is the thing that decides it.</p></figcaption>
    __C3__
  </figure>
</div>
<p class="col mu" style="margin-top:14px;font-size:.85rem">Each dot is a trait; vertical axis is how
strongly the correct tissue scored. Left: tracks available for that tissue (no correlation,
&rho;&nbsp;=&nbsp;&minus;0.25). Right: whether the model predicted any effect for that trait's
variants, measured without reference to any particular tissue (&rho;&nbsp;=&nbsp;+0.87).</p>

<section class="col stack">
  <h2>The clinching comparison</h2>
  <p class="callout">The <strong>brain</strong> is one of the best-measured tissues in the model
  &mdash; 63 datasets. Both brain traits fail.</p>
  <p class="callout good">The <strong>pancreatic islet</strong> is the worst-measured &mdash; 5
  datasets, and no open-chromatin data at all. Diabetes still ranks it first.</p>
  <p>So the question isn't "has the model seen this tissue?" It's <strong>"can the model see
  anything happening at these variants?"</strong> When it can, the tissue it names is usually
  right. When it can't, the ranking is noise &mdash; and nothing in the output tells you which
  situation you're in. That's the practical warning.</p>
</section>

<section class="col stack">
  <h2>So what about diabetes?</h2>
  <p>Honestly: <strong>this method can't settle it.</strong> In the larger study the islet ranks
  first, which matches decades of pancreas research. In a second, independent study it ranks 8th.
  The two studies' variant lists overlap by only 15 variants out of ~130 &mdash; they're nearly
  different samples of the same disease.</p>
  <p>Diabetes sits in the low-visibility zone where this benchmark is uninformative in either
  direction. Not evidence for the pancreas, not evidence against it.</p>
</section>

<section class="col">
  <div class="scroll">
  <table>
    <caption>Full results. Rank is the position of the known answer among 11 tissues; q is the
    false-discovery-corrected significance of that tissue against background.</caption>
    <thead><tr><th>Trait</th><th>Known tissue</th><th class="n">Rank</th><th class="n">Score</th><th class="n">q</th></tr></thead>
    <tbody>__TAB__</tbody>
  </table>
  </div>
</section>

<section class="col stack">
  <h2>What I'd do next</h2>
  <ul class="plain">
    <li>Compare against the standard statistical method for this question, to see whether
      per-variant prediction adds anything at all.</li>
    <li>Use whole credible sets rather than one variant per signal &mdash; the diabetes instability
      suggests a lot of the noise is upstream, in the fine-mapping.</li>
    <li>Test the "visibility" rule going forward: predict which traits will work <em>before</em>
      looking at the tissue answers.</li>
  </ul>
</section>

<footer class="col">
  <p>Predictions, not measurements. Nothing here was tested at a bench. Full method, statistics and
  citations are in the accompanying paper.</p>
</footer>
</div>
<div id="tip" role="status"></div>
"""
script = open(f"{SP}/tpl.html").read()
script = script[script.index("<script>"):]
html = ("<title>When Is a Tissue Prediction Trustworthy?</title>\n"+P['fonts']+"\n"+P['style']+"\n"
        + body.replace("__C1__",P['c1']).replace("__C2__",P['c2']).replace("__C3__",P['c3']).replace("__TAB__",P['tab'])
        + script)
open(f"{SP}/report.html","w").write(html)
print("wrote",len(html),"bytes")
