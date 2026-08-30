import { useState } from "react";
// import { LiteratureSearchBox } from "./LiteratureSearchBox"; // Literature tab disabled -- agent turn budget too tight, see rag/literature_agent.py.
// import { MethodologyQABox } from "./MethodologyQABox"; // Methodology tab disabled -- not evaluated/stabilized yet.
import { NarratorPanel } from "./NarratorPanel";
import type { RankResponse } from "../types/scoring";

type AssistantTab = "narrator" | "methodology" | "literature";

interface AiAssistantWidgetProps {
  result: RankResponse | null;
}

// Floating, minimized-by-default home for the read-only AI features (result narration;
// methodology Q&A and literature search are currently disabled, see tabs below) --
// QueryExpansionPanel stays embedded in QueryForm since it directly edits the inclusion gene list
// rather than producing read-only output. The panel and its tabs stay mounted for the widget's
// whole lifetime (only CSS-hidden when closed or inactive) so minimizing the widget or switching
// tabs never discards in-progress or completed question text/answers.
export function AiAssistantWidget({ result }: AiAssistantWidgetProps) {
  const [isOpen, setIsOpen] = useState(false);
  const [activeTab, setActiveTab] = useState<AssistantTab>("narrator");

  return (
    <div className="ai-assistant">
      <div className={isOpen ? "ai-assistant-panel" : "ai-assistant-panel ai-assistant-panel-hidden"}>
        <div className="ai-assistant-tabs">
          <button
            type="button"
            className={activeTab === "narrator" ? "active" : ""}
            onClick={() => setActiveTab("narrator")}
          >
            Narrator
          </button>
          {/* Methodology tab disabled -- not evaluated/stabilized yet.
          <button
            type="button"
            className={activeTab === "methodology" ? "active" : ""}
            onClick={() => setActiveTab("methodology")}
          >
            Methodology
          </button>
          */}
          {/* Literature tab disabled -- agent turn budget too tight, see rag/literature_agent.py.
          <button
            type="button"
            className={activeTab === "literature" ? "active" : ""}
            onClick={() => setActiveTab("literature")}
          >
            Literature
          </button>
          */}
        </div>
        <div className="ai-assistant-body">
          <div className={activeTab === "narrator" ? "" : "ai-assistant-tab-hidden"}>
            {result ? (
              <NarratorPanel evidenceRecord={result} />
            ) : (
              <p className="muted">Run a ranking query to enable the AI narrative.</p>
            )}
          </div>
          {/* Methodology tab disabled -- not evaluated/stabilized yet.
          <div className={activeTab === "methodology" ? "" : "ai-assistant-tab-hidden"}>
            <MethodologyQABox />
          </div>
          */}
          {/* Literature tab disabled -- agent turn budget too tight, see rag/literature_agent.py.
          <div className={activeTab === "literature" ? "" : "ai-assistant-tab-hidden"}>
            <LiteratureSearchBox />
          </div>
          */}
        </div>
      </div>
      <button
        type="button"
        className="ai-assistant-fab"
        aria-label={isOpen ? "Close AI assistant" : "Open AI assistant"}
        onClick={() => setIsOpen((open) => !open)}
      >
        {isOpen ? "×" : "AI"}
      </button>
    </div>
  );
}
