const copyButtonIcon = `
  <svg aria-hidden="true" viewBox="0 0 24 24" focusable="false">
    <path d="M16 1H6C4.9 1 4 1.9 4 3V17H6V3H16V1ZM19 5H10C8.9 5 8 5.9 8 7V21C8 22.1 8.9 23 10 23H19C20.1 23 21 22.1 21 21V7C21 5.9 20.1 5 19 5ZM19 21H10V7H19V21Z"></path>
  </svg>
`;

const trimCodeBlocks = () => {
  Array.from(document.getElementsByTagName("code")).forEach((block) => {
    block.innerHTML = block.innerHTML.trim();
  });
};

const getCopyText = (button) => {
  const copyTarget = button.dataset.copyTarget;

  if (copyTarget) {
    const target = document.querySelector(copyTarget);

    if (target) {
      return (target.innerText || target.textContent || "").trim();
    }
  }

  const codeBlock = button.closest("pre")?.querySelector("code");

  if (codeBlock) {
    return (codeBlock.innerText || codeBlock.textContent || "").trim();
  }

  return "";
};

const setCopyButtonState = (button, state = "default") => {
  const defaultLabel = button.dataset.copyLabel || "Copy";
  const successLabel = button.dataset.copySuccessLabel || "Copied";
  const errorLabel = button.dataset.copyErrorLabel || "Copy failed";
  const label =
    state === "copied"
      ? successLabel
      : state === "error"
        ? errorLabel
        : defaultLabel;

  button.classList.toggle("copied", state === "copied");
  button.classList.toggle("copy-error", state === "error");
  button.setAttribute("aria-label", label);
  button.setAttribute("title", label);

  const text = button.querySelector(".copy-button-label");

  if (text) {
    text.textContent = label;
  }
};

const fallbackCopyText = (text) => {
  const textArea = document.createElement("textarea");
  textArea.value = text;
  textArea.setAttribute("readonly", "");
  textArea.className = "sr-only";
  document.body.appendChild(textArea);
  textArea.select();
  const copied = document.execCommand("copy");
  document.body.removeChild(textArea);

  if (!copied) {
    throw new Error("Copy command failed");
  }
};

const handleCopyClick = async (event) => {
  const button = event.currentTarget;
  const text = getCopyText(button);

  if (!text) {
    return;
  }

  try {
    if (navigator.clipboard && window.isSecureContext) {
      await navigator.clipboard.writeText(text);
    } else {
      fallbackCopyText(text);
    }

    setCopyButtonState(button, "copied");
    window.setTimeout(() => setCopyButtonState(button), 2000);
  } catch {
    setCopyButtonState(button, "error");
    window.setTimeout(() => setCopyButtonState(button), 2000);
  }
};

const initializeCopyButton = (button) => {
  if (button.dataset.copyButtonReady === "true") {
    return;
  }

  button.type = "button";
  button.classList.add("copy-button");

  if (!button.querySelector(".copy-button-label")) {
    button.innerHTML = `${copyButtonIcon}<span class="copy-button-label">${button.dataset.copyLabel || "Copy"}</span>`;
  }

  setCopyButtonState(button);
  button.addEventListener("click", handleCopyClick);
  button.dataset.copyButtonReady = "true";
};

const addCodeCopyButton = (pre) => {
  if (pre.querySelector(".copy-button")) {
    return;
  }

  const button = document.createElement("button");
  button.className = "copy-button code-copy-button";
  button.dataset.copyLabel = "Copy";
  button.dataset.copySuccessLabel = "Copied";
  pre.appendChild(button);
  initializeCopyButton(button);
};

const initializeCopyButtons = () => {
  document.querySelectorAll("pre[data-copyable]").forEach(addCodeCopyButton);
  document.querySelectorAll(".copy-button[data-copy-target]").forEach(initializeCopyButton);
};

document.addEventListener("DOMContentLoaded", () => {
  trimCodeBlocks();
  initializeCopyButtons();
  window.hljs.highlightAll();
});

window.initializeCopyButtons = initializeCopyButtons;
