import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { Search, Database, Download, Upload, BarChart2, Dna, Beaker } from 'lucide-react';

const HomePage = () => {
  const [searchQuery, setSearchQuery] = useState('');
  const navigate = useNavigate();

  const handleSearch = (e: React.FormEvent) => {
    e.preventDefault();
    if (searchQuery.trim()) {
      navigate(`/search?q=${encodeURIComponent(searchQuery.trim())}`);
    }
  };

  return (
    <div className="min-h-screen bg-amber-50 flex flex-col">
      <div className="flex-grow">
        {/* Hero Section */}
        <section className="bg-gradient-to-r from-amber-500 to-red-500 text-white py-16">
          <div className="max-w-7xl mx-auto px-4">
            <div className="flex items-center justify-center mb-6">
              <div className="relative">
                <Beaker className="w-16 h-16 text-white" />
                <Dna className="w-8 h-8 text-white absolute -bottom-1 -right-1" />
              </div>
            </div>
            <h1 className="text-4xl font-bold mb-4 text-center">HYPLPDB</h1>
            <p className="text-xl text-center">Hypolipidemic Peptide Database</p>
          </div>
        </section>

        {/* Information Section */}
        <section className="max-w-7xl mx-auto px-4 py-12">
          <div className="space-y-12">
            {/* Content pairs */}
            <div className="flex flex-col md:flex-row gap-8 items-center">
              <div className="w-full md:w-1/2">
                <div className="bg-white rounded-2xl p-8 shadow-lg h-full">
                  <h2 className="text-2xl font-bold text-orange-600 mb-4">Understanding High Blood Lipids</h2>
                  <p className="text-gray-700 leading-relaxed">
                    High blood lipids can lead to serious health complications including heart disease and stroke. The accumulation of excess lipids in blood vessels creates blockages that impair circulation and increase health risks.
                  </p>
                </div>
              </div>
              <div className="w-full md:w-1/2">
                <img 
                  src="https://images.unsplash.com/photo-1630384060421-cb20d0e0649d?auto=format&fit=crop&q=80&w=800"
                  alt="Fried chicken" 
                  className="rounded-2xl shadow-lg object-cover w-full h-[250px]"
                />
              </div>
            </div>

            <div className="flex flex-col md:flex-row gap-8 items-center">
              <div className="w-full md:w-1/2">
                <div className="bg-white rounded-2xl p-8 shadow-lg h-full">
                  <h2 className="text-2xl font-bold text-orange-600 mb-4">Benefits of Lipid Management</h2>
                  <p className="text-gray-700 leading-relaxed">
                    Effective lipid management is crucial for cardiovascular health. Regular monitoring and control of blood lipid levels, combined with a healthy lifestyle, significantly reduces the risk of heart disease.
                  </p>
                </div>
              </div>
              <div className="w-full md:w-1/2">
                <img 
                  src="https://images.unsplash.com/photo-1606787366850-de6330128bfc?auto=format&fit=crop&q=80&w=800"
                  alt="French fries" 
                  className="rounded-2xl shadow-lg object-cover w-full h-[250px]"
                />
              </div>
            </div>

            <div className="flex flex-col md:flex-row gap-8 items-center">
              <div className="w-full md:w-1/2">
                <div className="bg-white rounded-2xl p-8 shadow-lg h-full">
                  <h2 className="text-2xl font-bold text-orange-600 mb-4">Innovative Peptide Solutions</h2>
                  <p className="text-gray-700 leading-relaxed">
                    Hypolipidemic peptides represent a breakthrough in lipid management, offering efficient and safe treatment options with minimal side effects. Ongoing research continues to unlock new potential in peptide-based therapies.
                  </p>
                </div>
              </div>
              <div className="w-full md:w-1/2">
                <img 
                  src="https://images.unsplash.com/photo-1529692236671-f1f6cf9683ba?auto=format&fit=crop&q=80&w=800"
                  alt="Fast food" 
                  className="rounded-2xl shadow-lg object-cover w-full h-[250px]"
                />
              </div>
            </div>
          </div>
        </section>

        {/* Search Section */}
        <section className="max-w-7xl mx-auto px-4 py-12 bg-white rounded-2xl shadow-lg mb-12">
          <h2 className="text-2xl font-bold text-orange-600 mb-6">Quick Search</h2>
          <form onSubmit={handleSearch} className="relative">
            <input
              type="text"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              placeholder="Search peptides..."
              className="w-full px-6 py-3 rounded-full border-2 border-amber-200 focus:border-amber-500 focus:ring-0"
            />
            <button
              type="submit"
              className="absolute right-2 top-1/2 -translate-y-1/2 px-6 py-2 bg-amber-500 text-white rounded-full hover:bg-amber-600 transition-colors"
            >
              <Search className="w-4 h-4" />
            </button>
          </form>
        </section>

        {/* Quick Links */}
        <section className="max-w-7xl mx-auto px-4 py-12">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
            <div 
              onClick={() => navigate('/analytics')}
              className="bg-white p-6 rounded-2xl shadow-lg hover:shadow-xl transition-shadow cursor-pointer"
            >
              <BarChart2 className="h-8 w-8 text-amber-500 mb-4" />
              <h3 className="text-xl font-bold text-gray-900 mb-2">Statistical Analysis</h3>
              <p className="text-gray-600">View statistical information and analysis reports of peptide data</p>
            </div>
            <div 
              onClick={() => navigate('/download')}
              className="bg-white p-6 rounded-2xl shadow-lg hover:shadow-xl transition-shadow cursor-pointer"
            >
              <Download className="h-8 w-8 text-amber-500 mb-4" />
              <h3 className="text-xl font-bold text-gray-900 mb-2">Data Download</h3>
              <p className="text-gray-600">Download complete datasets and structural information</p>
            </div>
            <div 
              onClick={() => navigate('/submit')}
              className="bg-white p-6 rounded-2xl shadow-lg hover:shadow-xl transition-shadow cursor-pointer"
            >
              <Upload className="h-8 w-8 text-amber-500 mb-4" />
              <h3 className="text-xl font-bold text-gray-900 mb-2">Submit Data</h3>
              <p className="text-gray-600">Submit new hypolipidemic peptide-related data</p>
            </div>
          </div>
        </section>
      </div>

      {/* Footer */}
      <footer className="bg-white py-6 mt-12">
        <div className="max-w-7xl mx-auto px-4">
          <div className="border-t border-gray-200 pt-6">
            <div className="text-center text-gray-600 text-sm">
              <p className="font-medium mb-2">HYPLPDB - Hypolipidemic Peptide Database</p>
              <p>© 2024 All rights reserved</p>
            </div>
          </div>
        </div>
      </footer>
    </div>
  );
};

export default HomePage;