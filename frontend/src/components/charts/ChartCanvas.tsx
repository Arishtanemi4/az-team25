import { useEffect, useRef } from "react";
import Chart, { type ChartConfiguration } from "chart.js/auto";

interface ChartCanvasProps {
  config: ChartConfiguration;
  height?: number;
}

// The one Chart.js wrapper in the app (first use of the library anywhere in the repo, V7-4
// pinned it, nothing had imported it yet). "chart.js/auto" auto-registers every
// controller/element/scale so call sites never need Chart.register(...) bookkeeping. Rebuilds
// the Chart instance whenever `config` changes (JSON-stringified so a fresh object literal each
// render doesn't retrigger it) and always destroys the previous instance first -- Chart.js
// leaks a canvas-bound instance otherwise.
export function ChartCanvas({ config, height = 320 }: ChartCanvasProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const chartRef = useRef<Chart | null>(null);
  const configKey = JSON.stringify(config);

  useEffect(() => {
    if (!canvasRef.current) return;
    chartRef.current?.destroy();
    chartRef.current = new Chart(canvasRef.current, config);
    return () => {
      chartRef.current?.destroy();
      chartRef.current = null;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [configKey]);

  return <canvas ref={canvasRef} height={height} />;
}
