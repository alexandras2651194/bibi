const TREE_URL = "../results/tree.newick";

const container = document.getElementById("treeContainer");
const status = document.getElementById("status");
const fitButton = document.getElementById("fitButton");

let currentTree = null;

function parseNewick(text) {
  const tokens = text.trim().replace(/;$/, "").split(/([(),:])/).map(t => t.trim()).filter(Boolean);
  let i = 0;

  function parseNode() {
    if (tokens[i] === "(") {
      i++;
      const children = [];
      while (true) {
        children.push(parseNode());
        if (tokens[i] === ",") { i++; continue; }
        if (tokens[i] === ")") { i++; break; }
        throw new Error("Invalid Newick structure near token " + tokens[i]);
      }

      let name = "";
      if (tokens[i] && ![",", ")", ":"].includes(tokens[i])) name = tokens[i++];
      let length = 0;
      if (tokens[i] === ":") { i++; length = Number(tokens[i++]); }
      return { name, length, children };
    }

    const name = tokens[i++];
    let length = 0;
    if (tokens[i] === ":") { i++; length = Number(tokens[i++]); }
    return { name, length, children: [] };
  }

  return parseNode();
}

function layoutTree(root) {
  const leaves = [];
  const internals = [];

  function visit(node, depth) {
    node.depth = depth;
    if (!node.children.length) {
      leaves.push(node);
      return;
    }
    internals.push(node);
    node.children.forEach(child => visit(child, depth + 1));
  }
  visit(root, 0);

  let leafIndex = 0;
  function assignY(node) {
    if (!node.children.length) {
      node.y = leafIndex++;
      return node.y;
    }
    node.children.forEach(assignY);
    node.y = node.children.reduce((sum, c) => sum + c.y, 0) / node.children.length;
    return node.y;
  }
  assignY(root);

  // A readable cladogram. Branch lengths are retained in the Newick source,
  // while the first visual version focuses on the tree topology/contours.
  const width = Math.max(1050, 720 + Math.max(0, root.depth) * 15);
  const rowHeight = 58;
  const top = 48;
  const height = Math.max(520, top * 2 + (leaves.length - 1) * rowHeight);

  function maxDepth(node) {
    if (!node.children.length) return node.depth;
    return Math.max(...node.children.map(maxDepth));
  }
  const depthMax = maxDepth(root) || 1;

  function xFor(node) {
    return 70 + (node.depth / depthMax) * (width - 350);
  }
  function yFor(node) {
    return top + node.y * rowHeight;
  }

  return { root, leaves, internals, width, height, xFor, yFor };
}

function esc(text) {
  return text.replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;").replaceAll('"', "&quot;");
}

function drawTree(root) {
  const layout = layoutTree(root);
  const { width, height, xFor, yFor } = layout;

  let svg = `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Neighbour Joining phylogenetic tree">`;
  svg += `<rect width="100%" height="100%" fill="#fbfcfd"/>`;

  function draw(node) {
    const x = xFor(node);
    const y = yFor(node);

    if (node.children.length) {
      const childYs = node.children.map(c => yFor(c));
      svg += `<line class="branch" x1="${x}" y1="${Math.min(...childYs)}" x2="${x}" y2="${Math.max(...childYs)}"/>`;
      node.children.forEach(child => {
        const cx = xFor(child);
        const cy = yFor(child);
        svg += `<line class="branch" x1="${x}" y1="${cy}" x2="${cx}" y2="${cy}"/>`;
        draw(child);
      });
      svg += `<circle class="node" cx="${x}" cy="${y}" r="4"/>`;
    } else {
      svg += `<circle class="node" cx="${x}" cy="${y}" r="4"/>`;
      svg += `<text class="leaf-label" x="${x + 12}" y="${y}">${esc(node.name)}</text>`;
    }
  }

  draw(root);
  svg += `</svg>`;
  container.innerHTML = svg;
  currentTree = root;
}

async function loadTree() {
  try {
    const response = await fetch(TREE_URL, { cache: "no-store" });
    if (!response.ok) throw new Error(`Could not load ${TREE_URL} (${response.status})`);
    const text = await response.text();
    const tree = parseNewick(text);
    drawTree(tree);
    const leaves = [];
    (function count(node) { if (!node.children.length) leaves.push(node); else node.children.forEach(count); })(tree);
    status.textContent = `${leaves.length} sequences · Newick loaded`;
  } catch (error) {
    container.innerHTML = `<div class="error"><strong>Tree unavailable.</strong><br>${esc(error.message)}<br><br>Run the C++ program and choose <code>T</code>, then reload this page through a local web server.</div>`;
    status.textContent = "Tree unavailable";
  }
}

fitButton.addEventListener("click", () => {
  if (currentTree) drawTree(currentTree);
});

loadTree();
