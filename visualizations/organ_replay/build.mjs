import { build, transform } from 'esbuild';
import { readFile, writeFile } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const here = path.dirname(fileURLToPath(import.meta.url));
const result = await build({ entryPoints: [path.join(here, 'viewer.js')], bundle: true, minify: true, format: 'iife', target: ['es2020'], write: false });
const template = await readFile(path.join(here, 'template.html'), 'utf8');
const styles = await readFile(path.join(here, 'style.css'), 'utf8');
const data = await readFile(path.join(here, 'data.json'), 'utf8');
const liver = (await readFile(path.join(here, 'assets/VH_M_Liver.glb'))).toString('base64');
const script = result.outputFiles[0].text.replaceAll('</script', '<\\/script');
const html = template.replace('/*__CABOMB_STYLE__*/', () => styles).replace('/*__CABOMB_DATA__*/', () => `window.CABOMB_DATA=${data};window.CABOMB_LIVER="${liver}";`).replace('/*__CABOMB_SCRIPT__*/', () => script);
for (const match of html.matchAll(/<script>([\s\S]*?)<\/script>/g)) {
  await transform(match[1], { loader: 'js', target: 'es2020' });
}
await writeFile(path.join(here, 'index.html'), html);
console.log(`Built offline index.html: ${(Buffer.byteLength(html) / 1024 / 1024).toFixed(2)} MiB`);
