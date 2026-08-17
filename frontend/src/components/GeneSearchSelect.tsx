import { useEffect, useRef, useState } from "react";
import { searchGenes } from "../api/geneServiceClient";
import type { GeneMatch } from "../types/scoring";

interface GeneSearchSelectProps {
  label: string;
  selected: GeneMatch | null;
  onChange: (gene: GeneMatch | null) => void;
}

// Single-gene version of GeneCombobox.tsx -- same debounced /genes/search call and the same
// browse-alphabetically-on-focus behaviour, but holds one gene instead of a chip list. Used by
// the relationships widgets (gene<->cell-line, gene<->disease), which each need exactly one
// gene, not the inclusion/exclusion multi-select GeneCombobox is built for.
export function GeneSearchSelect({ label, selected, onChange }: GeneSearchSelectProps) {
  const [query, setQuery] = useState("");
  const [matches, setMatches] = useState<GeneMatch[]>([]);
  const [isOpen, setIsOpen] = useState(false);
  const [highlightedIndex, setHighlightedIndex] = useState(0);
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (query.trim().length < 2) {
      setMatches([]);
      return;
    }
    const timer = setTimeout(() => {
      searchGenes(query, 20)
        .then((results) => {
          setMatches(results);
          setIsOpen(true);
          setHighlightedIndex(0);
        })
        .catch(() => setMatches([]));
    }, 300);
    return () => clearTimeout(timer);
  }, [query]);

  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (containerRef.current && !containerRef.current.contains(event.target as Node)) {
        setIsOpen(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  function selectGene(gene: GeneMatch) {
    onChange(gene);
    setQuery("");
    setMatches([]);
    setIsOpen(false);
  }

  function browseAlphabetically() {
    searchGenes("", 20)
      .then((results) => {
        setMatches(results);
        setIsOpen(true);
        setHighlightedIndex(0);
      })
      .catch(() => setMatches([]));
  }

  function handleFocus() {
    if (query.trim().length === 0) {
      browseAlphabetically();
      return;
    }
    if (matches.length > 0) setIsOpen(true);
  }

  function handleInputClick() {
    if (query.trim().length === 0 && !isOpen) browseAlphabetically();
  }

  function handleKeyDown(event: React.KeyboardEvent) {
    if (!isOpen || matches.length === 0) return;
    if (event.key === "ArrowDown") {
      event.preventDefault();
      setHighlightedIndex((i) => Math.min(i + 1, matches.length - 1));
    } else if (event.key === "ArrowUp") {
      event.preventDefault();
      setHighlightedIndex((i) => Math.max(i - 1, 0));
    } else if (event.key === "Enter") {
      event.preventDefault();
      selectGene(matches[highlightedIndex]);
    } else if (event.key === "Escape") {
      setIsOpen(false);
    }
  }

  return (
    <div className="combobox" ref={containerRef}>
      <label className="combobox-label">{label}</label>
      {selected && (
        <div className="chip-row">
          <span className="chip">
            {selected.symbol}
            <button type="button" aria-label={`Clear ${selected.symbol}`} onClick={() => onChange(null)}>
              ×
            </button>
          </span>
        </div>
      )}
      <div className="combobox-input-row">
        <input
          type="text"
          value={query}
          placeholder="Search gene symbol or Ensembl ID..."
          onChange={(e) => setQuery(e.target.value)}
          onFocus={handleFocus}
          onClick={handleInputClick}
          onKeyDown={handleKeyDown}
        />
        <button
          type="button"
          className="combobox-caret"
          aria-label={`Browse ${label} alphabetically`}
          onClick={browseAlphabetically}
        >
          ▾
        </button>
      </div>
      {isOpen && query.trim().length === 0 && <p className="combobox-hint">browsing alphabetically — type to search</p>}
      {isOpen && matches.length > 0 && (
        <ul className="combobox-dropdown">
          {matches.map((gene, index) => (
            <li
              key={gene.ensembl_id}
              className={index === highlightedIndex ? "highlighted" : ""}
              onMouseEnter={() => setHighlightedIndex(index)}
              onClick={() => selectGene(gene)}
            >
              <strong>{gene.symbol}</strong> <span className="muted">{gene.ensembl_id}</span>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
