import { useEffect, useRef, useState } from "react";
import { searchGenes } from "../api/geneServiceClient";
import type { GeneMatch } from "../types/scoring";

interface GeneComboboxProps {
  label: string;
  selected: GeneMatch[];
  onChange: (genes: GeneMatch[]) => void;
  excluded: GeneMatch[]; // genes already picked in the *other* list (inclusion vs exclusion),
                          // so the same gene can't be added to both at once.
}

// Search+dropdown gene input, reused for both the inclusion and exclusion gene lists. Typing
// debounces a call to gene_service's /genes/search, so the 53,820-row gene_reference table is
// never held client-side -- only the current query's matches are.
export function GeneCombobox({ label, selected, onChange, excluded }: GeneComboboxProps) {
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
          const selectedIds = new Set(selected.map((g) => g.ensembl_id));
          const excludedIds = new Set(excluded.map((g) => g.ensembl_id));
          setMatches(results.filter((g) => !selectedIds.has(g.ensembl_id) && !excludedIds.has(g.ensembl_id)));
          setIsOpen(true);
          setHighlightedIndex(0);
        })
        .catch(() => setMatches([]));
    }, 300);
    return () => clearTimeout(timer);
  }, [query, selected, excluded]);

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
    onChange([...selected, gene]);
    setQuery("");
    setMatches([]);
    setIsOpen(false);
  }

  function removeGene(ensemblId: string) {
    onChange(selected.filter((g) => g.ensembl_id !== ensemblId));
  }

  // Shows the default alphabetical gene list (rather than nothing) so a researcher can browse
  // without already knowing a symbol. Shared by focus, click, and the caret button (V7-1d) --
  // clicking an already-focused input used to do nothing new, so the caret needs its own handler
  // that doesn't depend on a focus event having just fired.
  function browseAlphabetically() {
    searchGenes("", 20)
      .then((results) => {
        const selectedIds = new Set(selected.map((g) => g.ensembl_id));
        const excludedIds = new Set(excluded.map((g) => g.ensembl_id));
        setMatches(results.filter((g) => !selectedIds.has(g.ensembl_id) && !excludedIds.has(g.ensembl_id)));
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
      <div className="chip-row">
        {selected.map((gene) => (
          <span key={gene.ensembl_id} className="chip">
            {gene.symbol}
            <button type="button" aria-label={`Remove ${gene.symbol}`} onClick={() => removeGene(gene.ensembl_id)}>
              ×
            </button>
          </span>
        ))}
      </div>
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
