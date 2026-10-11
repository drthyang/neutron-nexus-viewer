// The NeXus Viewer as an MCP server: every agent tool (index.js) as an MCP tool, and the method notes as
// a resource. One Workspace per server, so the volumes a client opens stay open for its session.
// main.js runs it on stdio.

import { readFileSync } from 'node:fs';

import { McpServer } from '@modelcontextprotocol/sdk/server/mcp.js';

import { INSTRUCTIONS, runTool, TOOLS, toJson, Workspace } from './index.js';

const pkg = JSON.parse(readFileSync(new URL('../package.json', import.meta.url), 'utf8'));
export const VERSION = pkg.version ?? '0.0.0';

export function createServer(workspace = new Workspace()) {
  const server = new McpServer({ name: 'nexus-viewer', title: 'NeXus Viewer', version: VERSION }, { instructions: INSTRUCTIONS });
  for (const tool of TOOLS) {
    server.registerTool(tool.name, {
      title: tool.title, description: tool.description, inputSchema: tool.input, annotations: { title: tool.title, ...tool.annotations },
    }, async (args, extra) => {
      const r = await runTool(tool, args, { workspace, signal: extra.signal });
      if (!r.ok) return { content: [{ type: 'text', text: r.details ? `${r.error}\n${toJson(r.details)}` : r.error }], isError: true };
      const images = r.images.map((img) => ({ type: 'image', data: Buffer.from(img.data).toString('base64'), mimeType: img.mimeType }));
      return { content: [{ type: 'text', text: toJson(r.value) }, ...images] };
    });
  }
  server.registerResource(
    'method',
    'nexus-viewer://docs/method',
    { title: 'NeXus Viewer method', description: 'How slices, line cuts, symmetry averaging, masks, I(Q), the NEBULA3D export and the Rigaku reduction are computed.', mimeType: 'text/markdown' },
    (uri) => ({ contents: [{ uri: uri.href, mimeType: 'text/markdown', text: readFileSync(new URL('../docs/METHOD.md', import.meta.url), 'utf8') }] }),
  );
  return server;
}
