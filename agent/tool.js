// The NeXus Viewer's analysis as agent tools: a name, a description written for a model, a zod input
// schema and a handler that returns plain JSON. Nothing here depends on a transport: server.js serves
// the tools over MCP, and toolDefinitions() (index.js) gives them as JSON Schema for any tool-use API.

/** A tool, as { name, title, description, input (zod object), annotations, run(args, ctx) }. */
export const defineTool = (tool) => tool;

/** A failure the caller can act on: the message says what to change, `details` carries the choices (ids, presets). */
export class ToolError extends Error {
  constructor(message, details = undefined) {
    super(message);
    this.name = 'ToolError';
    this.details = details;
  }
}

// Images a result carries besides its JSON (PNG bytes), under a key JSON.stringify skips.
const IMAGES = Symbol('images');

/** `value` with images attached: [{ data: Uint8Array, mimeType }]. */
export function withImages(value, images) {
  Object.defineProperty(value, IMAGES, { value: images, enumerable: false });
  return value;
}

/** The images attached to a result, or []. */
export const imagesOf = (value) => value?.[IMAGES] ?? [];

/**
 * JSON for a tool result: numbers to 7 significant figures (integers unchanged), magnitudes below
 * 1e-12 as 0, typed arrays as arrays, and non-finite numbers as null.
 */
export function toJson(value) {
  return JSON.stringify(value, (_key, v) => {
    if (typeof v === 'number') return Number.isFinite(v) ? (Number.isInteger(v) ? v : Math.abs(v) < 1e-12 ? 0 : Number(v.toPrecision(7))) : null;
    if (ArrayBuffer.isView(v) && !(v instanceof DataView)) return Array.from(v);
    return v;
  });
}
