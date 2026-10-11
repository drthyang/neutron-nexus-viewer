#!/usr/bin/env node
// Runs the NeXus Viewer MCP server on stdio (`npm run mcp`). Logs go to stderr: stdout carries the protocol.
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';

import { createServer, VERSION } from './server.js';

await createServer().connect(new StdioServerTransport());
console.error(`NeXus Viewer ${VERSION} MCP server on stdio`);
