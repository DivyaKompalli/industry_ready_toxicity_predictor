// JS for interactive charts, filtering, and note saving
// Handle lead optimization form
document
  .getElementById("optimize-form")
  .addEventListener("submit", async (e) => {
    e.preventDefault();

    const formData = new FormData(e.target);
    const response = await fetch("/optimize", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        smiles: formData.get("smiles"),
        properties: {
          toxicity: parseFloat(formData.get("max_toxicity")),
          logp: parseFloat(formData.get("target_logp")),
        },
      }),
    });

    const results = await response.json();
    renderOptimizationResults(results);
  });

function renderOptimizationResults(results) {
  const container = document.getElementById("optimization-results");
  container.innerHTML = "";

  results.forEach((mol) => {
    const card = document.createElement("div");
    card.className = "molecule-card";
    card.innerHTML = `
            <img src="data:image/png;base64,${mol.image}" />
            <div class="molecule-props">
                <span>Toxicity: ${mol.toxicity.toFixed(2)}</span>
                <span>LogP: ${mol.logp.toFixed(2)}</span>
            </div>
        `;
    container.appendChild(card);
  });
}

// Initialize all tooltips
document.querySelectorAll("[data-tooltip]").forEach((el) => {
  tippy(el, {
    content: el.getAttribute("data-tooltip"),
    placement: "top",
  });
});
