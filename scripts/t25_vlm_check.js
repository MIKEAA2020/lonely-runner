const ZAI = require('z-ai-web-dev-sdk').default;
const fs = require('fs');

async function main() {
  const zai = await ZAI.create();
  const pages = ['t25_p49.png', 't25_p50.png', 't25_p52.png', 't25_p53.png'];
  const content = [{
    type: 'text',
    text: 'Strict layout QA of these four pages from a compiled mathematics paper ' +
          '(they contain a new lemma, its proof, and prose analysis). For EACH page ' +
          'report: (1) text running past the margins or cut off; (2) overlapping ' +
          'text or blocks; (3) garbled or missing mathematical notation; ' +
          '(4) any rendering artifact. End with a per-page PASS/FAIL verdict.'
  }];
  for (const f of pages) {
    const b64 = fs.readFileSync('/home/z/my-project/scripts/paper/' + f).toString('base64');
    content.push({ type: 'image_url', image_url: { url: 'data:image/png;base64,' + b64 } });
  }
  const r = await zai.chat.completions.createVision({
    messages: [{ role: 'user', content }],
    model: 'glm-5v-turbo'
  });
  const out = r.choices[0].message.content;
  fs.writeFileSync('/home/z/my-project/scripts/paper/t25_vlm_check.json', JSON.stringify(r, null, 1));
  console.log(out);
}
main().catch(e => { console.error('VLM check failed:', e.message); process.exit(1); });
