import { NavLink } from "react-router-dom";

// Every route CellLineSelector exposes today, in the order they appear in the nav bar.
const ROUTES: Array<{ to: string; label: string }> = [
  { to: "/", label: "Rank" },
  { to: "/data/schema", label: "Data schema" },
  { to: "/data/query", label: "Data query" },
  { to: "/explore/relationships", label: "Relationships" },
  { to: "/explore/eda", label: "EDA" },
  { to: "/validation", label: "Validation" },
  { to: "/graph", label: "Knowledge graph" },
];

export function NavBar() {
  return (
    <nav className="top-nav">
      <ul className="top-nav-list">
        {ROUTES.map((route) => (
          <li key={route.to}>
            <NavLink
              to={route.to}
              end={route.to === "/"}
              className={({ isActive }) => (isActive ? "top-nav-link top-nav-link-active" : "top-nav-link")}
            >
              {route.label}
            </NavLink>
          </li>
        ))}
      </ul>
    </nav>
  );
}
