// Write docs/screenshot.png: the example (?demo: A and B in split view, 6/mmm, the edge mask)
// with a line cut along (0.5, 0, L) in the fourth view, at 1440×900 and 1.5× pixel ratio, with
// "© 2026 Tsung-Han Yang" in the bottom-left corner.
// It starts a headless Chrome with a temporary profile and drives it over the DevTools protocol.
//
//   python3 -m http.server 8000 &
//   node tools/make_screenshot.mjs http://localhost:8000/?demo docs/screenshot.png
//
// CHROME may point at the browser binary; the default is Google Chrome on macOS.
import { spawn } from 'node:child_process';
import { mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const [url = 'http://localhost:8000/?demo', out = 'docs/screenshot.png'] = process.argv.slice(2);
const chrome = process.env.CHROME ?? '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const port = 9333, profile = mkdtempSync(join(tmpdir(), 'nxv-shot-'));
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const browser = spawn(chrome, ['--headless=new', `--remote-debugging-port=${port}`, `--user-data-dir=${profile}`,
  '--no-first-run', '--no-default-browser-check', '--hide-scrollbars', 'about:blank'], { stdio: 'ignore' });

try {
  let target;
  for (let i = 0; i < 50 && !target; i++) {
    await sleep(200);
    target = await fetch(`http://127.0.0.1:${port}/json/list`).then((r) => r.json()).then((l) => l.find((t) => t.type === 'page')).catch(() => null);
  }
  if (!target) throw new Error('Chrome did not start; set CHROME to its binary.');
  const ws = new WebSocket(target.webSocketDebuggerUrl);
  await new Promise((r) => ws.addEventListener('open', r, { once: true }));
  let next = 0;
  const pending = new Map();
  ws.addEventListener('message', ({ data }) => {
    const msg = JSON.parse(data);
    if (pending.has(msg.id)) { pending.get(msg.id)(msg); pending.delete(msg.id); }
  });
  const send = (method, params = {}) => new Promise((resolve) => {
    pending.set(++next, resolve);
    ws.send(JSON.stringify({ id: next, method, params }));
  });
  const run = async (expression) => (await send('Runtime.evaluate', { expression, awaitPromise: true, returnByValue: true })).result?.result?.value;

  await send('Emulation.setDeviceMetricsOverride', { width: 1440, height: 900, deviceScaleFactor: 1.5, mobile: false });
  await send('Page.navigate', { url });
  // Wait for A and B, the symmetry and the mask.
  let ready = false;
  for (let i = 0; i < 120 && !ready; i++) {
    await sleep(500);
    ready = await run(`!!document.querySelector('.ds-b:not([hidden])') && document.getElementById('mask-state')?.textContent.includes('removed')
      && document.getElementById('status-text')?.textContent === 'ready'`);
  }
  if (!ready) throw new Error(`The example did not finish loading from ${url}.`);
  await run(`(async () => {
    [...document.querySelectorAll('.slot-switch button')].find((b) => b.dataset.value === 'cut').click();
    await new Promise((r) => setTimeout(r, 500));
    for (const [sel, v] of [['.cut-from', '0.5, 0, -4'], ['.cut-to', '0.5, 0, 4'], ['.cut-width', '0.52']]) {
      const el = document.querySelector(sel);
      el.value = v;
      el.dispatchEvent(new Event('change'));
    }
  })()`);
  await sleep(4000);
  // Every demo image carries the copyright in its bottom-left corner, as in the other apps' media. It
  // takes the place of the panel's own credit line (its text is hidden, its bar kept), which it would overlap.
  await run(`(() => {
    document.querySelectorAll('.panel-foot > *').forEach((el) => { el.style.visibility = 'hidden'; });
    const mark = document.createElement('div');
    mark.textContent = '© 2026 Tsung-Han Yang';
    mark.style.cssText = 'position:fixed;left:8px;bottom:8px;z-index:2147483647;padding:2px 8px;'
      + 'background:rgba(0,0,0,0.72);color:#fff;font:15px/1.4 system-ui,-apple-system,"Segoe UI",sans-serif;pointer-events:none';
    document.body.append(mark);
  })()`);
  await sleep(200);
  const shot = await send('Page.captureScreenshot', { format: 'png' });
  writeFileSync(out, Buffer.from(shot.result.data, 'base64'));
  console.log(`wrote ${out}`);
  ws.close();
} finally {
  browser.kill();
  await sleep(500);
  rmSync(profile, { recursive: true, force: true });
}
