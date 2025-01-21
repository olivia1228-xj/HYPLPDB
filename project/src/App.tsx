import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Layout from './components/Layout';
import HomePage from './pages/HomePage';
import SearchPage from './pages/SearchPage';
import AnalyticsPage from './pages/AnalyticsPage';
import DownloadPage from './pages/DownloadPage';
import SubmitPage from './pages/SubmitPage';
import PeptideInfoPage from './pages/PeptideInfoPage';

function App() {
  return (
    <Router>
      <Layout>
        <Routes>
          <Route path="/" element={<HomePage />} />
          <Route path="/search" element={<SearchPage />} />
          <Route path="/analytics" element={<AnalyticsPage />} />
          <Route path="/download" element={<DownloadPage />} />
          <Route path="/submit" element={<SubmitPage />} />
          <Route path="/peptide/:id" element={<PeptideInfoPage />} />
        </Routes>
      </Layout>
    </Router>
  );
}

export default App;