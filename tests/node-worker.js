// Web-Worker-like wrapper around node:worker_threads, so the tests can drive the module
// workers of the app (postMessage / onmessage / onerror / terminate, transfer lists).
import { Worker } from 'node:worker_threads';

const SHIM = `
import { parentPort, workerData } from 'node:worker_threads';
globalThis.self = globalThis;
self.postMessage = (m, t) => parentPort.postMessage(m, t);
parentPort.on('message', (data) => self.onmessage?.({ data }));
await import(workerData.url);
`;

export function nodeWorker(url) {
  const w = new Worker(new URL(`data:text/javascript,${encodeURIComponent(SHIM)}`), { workerData: { url: String(url) } });
  const o = {
    onmessage: null, onerror: null,
    postMessage: (m, t) => w.postMessage(m, t),
    terminate: () => w.terminate(),
  };
  w.on('message', (data) => o.onmessage?.({ data }));
  w.on('error', (e) => o.onerror?.(e));
  return o;
}
