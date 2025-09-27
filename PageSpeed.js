const psi = require("psi");
const process = require("process");
const parseInputs = (inputs) => {
  const keyValuePairs = inputs.split(",");

  let url, strategy, threshold, apiKey;

  keyValuePairs.forEach((pair) => {
    const [key, value] = pair.split("=");
    switch (key) {
      case "url":
        url = value;
        break;
      case "strategy":
        strategy = value;
        break;
      case "threshold":
        threshold = value;
        break;
      case "apiKey":
        apiKey = value;
        break;
      default:
        throw new Error(
          "Invalid input format, please use the following: url=<your url>,strategy=<desktop/mobile>,threshold=<number(1-100)/none>,apiKey=<yourkey/none"
        );
    }
  });
  if (!url || !strategy || !threshold || !apiKey) {
    throw new Error(
      "Invalid input format,please use the following: url=<your url> strategy=<desktop/mobile>,threshold=<number 1-100/none>,apiKey=<yourkey/none"
    );
  }

  if (strategy !== "desktop" && strategy !== "mobile") {
    throw new Error('Invalid strategy. Enter "desktop" or "mobile".');
  }

  if (threshold !== "none") {
    const numThreshold = Number(threshold);
    if (isNaN(numThreshold) || numThreshold < 0 || numThreshold > 100) {
      throw new Error(
        'Invalid threshold. Enter a number between 0 and 100 or "none".'
      );
    }
  }

  return {
    url,
    strategy,
    threshold: threshold === "none" ? 70 : Number(threshold),
    key: apiKey === "none" ? undefined : apiKey,
  };
};

const runPageSpeed = async (inputs) => {
  const { url, strategy, threshold, key } = { ...inputs };
  try {
    console.log(`Page Speed results for ${url} using ${strategy}`);
    await psi.output(url, {
      ...(key ? { key } : undefined),
      ...(key ? undefined : { nokey: "true" }),
      strategy,
      format: "cli",
      threshold,
    });
  } catch (e) {
    throw e.message;
  }
};

let parsedInputs;
try {
  parsedInputs = parseInputs(process.argv[2]);
} catch (e) {
  throw e.message;
}

runPageSpeed(parsedInputs);
