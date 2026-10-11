import assert from 'node:assert/strict';
import test from 'node:test';

import { Client } from '@modelcontextprotocol/sdk/client/index.js';
import { InMemoryTransport } from '@modelcontextprotocol/sdk/inMemory.js';

import { TOOLS } from '../agent/index.js';
import { createServer } from '../agent/server.js';

const demo = new URL('../examples/demo_10K.nxs', import.meta.url).pathname;
const text = (r) => r.content.find((c) => c.type === 'text').text;

async function connect() {
  const client = new Client({ name: 'test', version: '0' });
  const [a, b] = InMemoryTransport.createLinkedPair();
  await createServer().connect(a);
  await client.connect(b);
  return client;
}

test('the MCP server lists every tool, keeps volumes open, and returns images', async () => {
  const client = await connect();
  try {
    assert.equal(client.getServerVersion().name, 'nexus-viewer');
    assert.match(client.getInstructions(), /open_volume/);
    const { tools } = await client.listTools();
    assert.deepEqual(tools.map((t) => t.name), TOOLS.map((t) => t.name));
    const slice = tools.find((t) => t.name === 'slice_volume');
    assert.deepEqual(slice.inputSchema.required, ['volume_id', 'fixed_axis']);

    const opened = JSON.parse(text(await client.callTool({ name: 'open_volume', arguments: { path: demo, id: 'cold' } })));
    assert.equal(opened.volume_id, 'cold');
    // The demo's raw volume has no data at K = 0: 6/mmm fills it from equivalent voxels.
    await client.callTool({ name: 'set_symmetry', arguments: { volume_id: 'cold', laue_class: '6/mmm' } });
    const r = await client.callTool({ name: 'slice_volume', arguments: { volume_id: 'cold', fixed_axis: 'K', center: 0 } });
    assert.ok(!r.isError);
    assert.equal(JSON.parse(text(r)).plane, 'HL');
    const image = r.content.find((c) => c.type === 'image');
    assert.equal(image.mimeType, 'image/png');
    assert.ok(Buffer.from(image.data, 'base64').length > 1000);

    const bad = await client.callTool({ name: 'line_cut', arguments: { volume_id: 'warm', start: [0, 0, 0], end: [1, 0, 0] } });
    assert.equal(bad.isError, true);
    assert.match(text(bad), /No volume 'warm'/);

    const method = await client.readResource({ uri: 'nexus-viewer://docs/method' });
    assert.match(method.contents[0].text, /^# Method/);
  } finally {
    await client.close();
  }
});
