import { Route, Routes } from "react-router-dom";
import { NavBar } from "./components/NavBar";
import DataQueryPage from "./pages/DataQueryPage";
import DataSchemaPage from "./pages/DataSchemaPage";
import EdaPage from "./pages/EdaPage";
import GraphPage from "./pages/GraphPage";
import RankPage from "./pages/RankPage";
import RelationshipsPage from "./pages/RelationshipsPage";
import ValidationPage from "./pages/ValidationPage";
import "./App.css";

function App() {
  return (
    <>
      <NavBar />
      <Routes>
        <Route path="/" element={<RankPage />} />
        <Route path="/data/schema" element={<DataSchemaPage />} />
        <Route path="/data/query" element={<DataQueryPage />} />
        <Route path="/explore/relationships" element={<RelationshipsPage />} />
        <Route path="/explore/eda" element={<EdaPage />} />
        <Route path="/validation" element={<ValidationPage />} />
        <Route path="/graph" element={<GraphPage />} />
      </Routes>
    </>
  );
}

export default App;
