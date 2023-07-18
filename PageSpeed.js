const psi = require("psi");
const process = require("process");
const url = process.argv[2];
const strategyInput = process.argv[3];
const thresholdInput = process.argv[4];
const keyInput = process.argv[5];
const run = async () => {
  try {
    const strategy = strategyInput?.toLowerCase();
    const key = keyInput || undefined;
    const threshold = Number(thresholdInput);
    if (!url) {
      throw new Error("A valid Url is required to run Page Speed Insights.");
    }
    if (!strategy || (strategy !== "desktop" && strategy !== "mobile")) {
      throw new Error(
        "A valid strategy is required to run Page Speed Insights. (desktop or mobile)"
      );
    }

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

run();
