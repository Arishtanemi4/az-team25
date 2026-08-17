import type { RankStageEvent } from "../api/scoringServiceClient";

interface RankingProgressProps {
  stages: RankStageEvent[];
}

// Replaces the old static "this can take a few minutes" sentence with a live checklist driven
// by the /rank SSE stream's `stage` events, so a researcher can see what's actually happening
// during a multi-minute ranking run instead of staring at a blank line.
export function RankingProgress({ stages }: RankingProgressProps) {
  const total = stages.at(-1)?.total ?? 8;

  return (
    <div className="ranking-progress">
      <p className="loading-notice">
        Scoring candidates against the full evidence tables -- this can take a few minutes.
      </p>
      <ol className="ranking-progress-steps">
        {stages.map((stage, i) => (
          <li key={stage.index} className={i === stages.length - 1 ? "current" : "done"}>
            {stage.label}
          </li>
        ))}
        {stages.length < total && <li className="pending">Waiting for the next step...</li>}
      </ol>
    </div>
  );
}
