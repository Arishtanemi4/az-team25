import type { DiagnosticInfo } from "../types/scoring";

interface DiagnosticBannerProps {
  diagnostic: DiagnosticInfo;
}

// Shown when every candidate scored D=0 (an over-constrained query) -- an empty ranked list is
// explained, never silently blank (ALGORITHM_SPEC.md edge case 5).
export function DiagnosticBanner({ diagnostic }: DiagnosticBannerProps) {
  return (
    <div className="diagnostic-banner">
      <p>{diagnostic.message}</p>
      {diagnostic.most_disqualifying_gene && (
        <p>
          Most disqualifying gene: <strong>{diagnostic.most_disqualifying_gene[1]}</strong> (
          {diagnostic.most_disqualifying_gene[0]})
        </p>
      )}
    </div>
  );
}
