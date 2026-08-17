import { useEffect, useRef, useState } from "react";

interface FilterSelectProps {
  label: string;
  options: string[];
  value: string | null;
  onChange: (value: string | null) => void;
  placeholder?: string;
}

// Search+dropdown single-select for the lineage/primary_disease filters. Unlike GeneCombobox
// this filters a small, already-fetched list client-side (there are only dozens of lineages),
// so no debounce or server round-trip is needed.
export function FilterSelect({ label, options, value, onChange, placeholder }: FilterSelectProps) {
  const [query, setQuery] = useState("");
  const [isOpen, setIsOpen] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (containerRef.current && !containerRef.current.contains(event.target as Node)) {
        setIsOpen(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  const filteredOptions = options.filter((o) => o.toLowerCase().includes(query.toLowerCase()));

  return (
    <div className="combobox" ref={containerRef}>
      <label className="combobox-label">{label}</label>
      <div className="filter-input-row">
        <input
          type="text"
          value={value ?? query}
          placeholder={placeholder ?? `Any ${label.toLowerCase()}`}
          onFocus={() => {
            setQuery("");
            setIsOpen(true);
          }}
          onChange={(e) => {
            onChange(null);
            setQuery(e.target.value);
          }}
        />
        {value && (
          <button type="button" className="clear-button" aria-label={`Clear ${label}`} onClick={() => onChange(null)}>
            ×
          </button>
        )}
      </div>
      {isOpen && (
        <ul className="combobox-dropdown">
          {filteredOptions.length === 0 && <li className="muted">No matches</li>}
          {filteredOptions.map((option) => (
            <li
              key={option}
              onClick={() => {
                onChange(option);
                setQuery("");
                setIsOpen(false);
              }}
            >
              {option}
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
