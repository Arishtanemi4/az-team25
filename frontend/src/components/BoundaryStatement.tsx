interface BoundaryStatementProps {
  text: string;
}

// Renders the API's boundary_statement verbatim -- root _.md SS2's confidence-vs-proof
// distinction must be visible wherever results are shown, not just in documentation.
export function BoundaryStatement({ text }: BoundaryStatementProps) {
  return <p className="boundary-statement">{text}</p>;
}
