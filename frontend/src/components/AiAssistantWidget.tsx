import { useState } from "react";
import { LiteratureSearchBox } from "./LiteratureSearchBox";
import { MethodologyQABox } from "./MethodologyQABox";
import { NarratorPanel } from "./NarratorPanel";
import type { RankResponse } from "../types/scoring";

type AssistantTab = "narrator" | "methodology" | "literature";

interface AiAssistantWidgetProps {
  result: RankResponse | null;
}

// Floating, minimized-by-default home for the three read-only AI features (result narration,
// methodology Q&A, literature search) -- QueryExpansionPanel stays embedded in QueryForm since
// it directly edits the inclusion gene list rather than producing read-only output. The panel and
// all three tabs stay mounted for the widget's whole lifetime (only CSS-hidden when closed or
// inactive) so minimizing the widget or switching tabs never discards in-progress or completed
// question text/answers.
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
          <button
            type="button"
            className={activeTab === "methodology" ? "active" : ""}
            onClick={() => setActiveTab("methodology")}
          >
            Methodology
          </button>
          <button
            type="button"
            className={activeTab === "literature" ? "active" : ""}
            onClick={() => setActiveTab("literature")}
          >
            Literature
          </button>
        </div>
        <div className="ai-assistant-body">
          <div className={activeTab === "narrator" ? "" : "ai-assistant-tab-hidden"}>
            {result ? (
              <NarratorPanel evidenceRecord={result} />
            ) : (
              <p className="muted">Run a ranking query to enable the AI narrative.</p>
            )}
          </div>
          <div className={activeTab === "methodology" ? "" : "ai-assistant-tab-hidden"}>
            <MethodologyQABox />
          </div>
          <div className={activeTab === "literature" ? "" : "ai-assistant-tab-hidden"}>
            <LiteratureSearchBox />
          </div>
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
