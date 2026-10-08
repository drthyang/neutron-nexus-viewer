// A small pool of workers (rigaku-map-worker.js) for the per-frame work of the Rigaku
// reduction (rigaku-reduce.js). Workers are anything with postMessage / onmessage / onerror /
// terminate: Web Workers in the browser, a worker_threads wrapper in the tests.

export class WorkerPool {
  constructor(workers) {
    this.workers = workers;
    this.load = workers.map(() => 0);
    this.pending = new Map();
    this.next = 1;
    this.broken = null;
    workers.forEach((w, i) => {
      w.onmessage = ({ data }) => {
        const p = this.pending.get(data.id);
        if (!p) return;
        this.pending.delete(data.id);
        this.load[i]--;
        if (data.error) p.reject(new Error(data.error));
        else p.resolve(data);
      };
      w.onerror = (e) => {
        e?.preventDefault?.();
        this.broken = new Error(e?.message || 'A reduction worker failed.');
        for (const p of this.pending.values()) p.reject(this.broken);
        this.pending.clear();
      };
    });
  }

  get size() { return this.workers.length; }

  /** Send `msg` to worker `i`; resolves with its reply. */
  call(i, msg, transfer = []) {
    if (this.broken) return Promise.reject(this.broken);
    const id = this.next++;
    this.load[i]++;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.workers[i].postMessage({ ...msg, id }, transfer);
    });
  }

  /** Send `msg` to every worker; resolves with their replies. */
  all(msg) { return Promise.all(this.workers.map((_, i) => this.call(i, msg))); }

  /** Resolves when every worker answers, or rejects after `ms` milliseconds. */
  ready(ms = 20000) {
    let timer;
    const timeout = new Promise((_, reject) => { timer = setTimeout(() => reject(new Error('The reduction workers did not start.')), ms); });
    return Promise.race([this.all({ type: 'ping' }), timeout]).finally(() => clearTimeout(timer));
  }

  /**
   * Tasks 0..n-1: `start(i)` (may be async) returns { msg, transfer } for task i, which goes to
   * the least busy worker; `apply(i, reply)` is called for the replies strictly in task order,
   * with at most `depth` tasks per worker in flight.
   */
  async ordered(n, start, apply, depth = 2) {
    const cap = this.workers.length * depth, queue = [];
    const settle = async () => { const q = queue.shift(); await apply(q.i, await q.p); };
    for (let i = 0; i < n; i++) {
      if (queue.length >= cap) await settle();
      const { msg, transfer } = await start(i);
      const w = this.load.indexOf(Math.min(...this.load));
      const p = this.call(w, msg, transfer);
      p.catch(() => {}); // awaited in order by settle(); avoids an unhandled-rejection report meanwhile
      queue.push({ i, p });
    }
    while (queue.length) await settle();
  }

  terminate() { for (const w of this.workers) w.terminate(); }
}
