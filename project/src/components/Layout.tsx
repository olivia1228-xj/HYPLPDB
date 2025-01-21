import React from 'react';
import { Link } from 'react-router-dom';
import { Database, Search, BarChart2, Download, Upload } from 'lucide-react';

const Layout = ({ children }: { children: React.ReactNode }) => {
  return (
    <div className="min-h-screen bg-amber-50">
      <nav className="bg-white shadow">
        <div className="max-w-7xl mx-auto px-4">
          <div className="flex items-center justify-between h-16">
            <Link to="/" className="flex items-center">
              <Database className="h-8 w-8 text-amber-500" />
              <span className="ml-2 text-xl font-bold text-gray-900">HYPLPDB</span>
            </Link>
            <div className="flex space-x-4">
              <Link to="/search" className="flex items-center text-gray-700 hover:text-amber-500">
                <Search className="h-5 w-5 mr-1" />
                Search
              </Link>
              <Link to="/analytics" className="flex items-center text-gray-700 hover:text-amber-500">
                <BarChart2 className="h-5 w-5 mr-1" />
                Analytics
              </Link>
              <Link to="/download" className="flex items-center text-gray-700 hover:text-amber-500">
                <Download className="h-5 w-5 mr-1" />
                Download
              </Link>
              <Link to="/submit" className="flex items-center text-gray-700 hover:text-amber-500">
                <Upload className="h-5 w-5 mr-1" />
                Submit
              </Link>
            </div>
          </div>
        </div>
      </nav>
      <main className="max-w-7xl mx-auto px-4 py-8">
        {children}
      </main>
    </div>
  );
};

export default Layout;