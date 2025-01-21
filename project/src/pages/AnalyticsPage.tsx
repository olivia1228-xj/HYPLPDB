import React, { useState, useEffect } from 'react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, PieChart, Pie, Cell } from 'recharts';
import { Download, Loader2, AlertCircle } from 'lucide-react';
import * as XLSX from 'xlsx';
import { getAnalyticsData } from '../lib/api';

const COLORS = [
  '#3B82F6', // Blue
  '#10B981', // Green
  '#F59E0B', // Yellow
  '#EF4444', // Red
  '#8B5CF6', // Purple
  '#EC4899', // Pink
  '#14B8A6', // Teal
  '#F97316', // Orange
  '#6366F1', // Indigo
  '#84CC16'  // Lime
];

const CustomPieLabel = ({ cx, cy, midAngle, innerRadius, outerRadius, value, name, percent }) => {
  const RADIAN = Math.PI / 180;
  const radius = outerRadius * 1.2;
  const x = cx + radius * Math.cos(-midAngle * RADIAN);
  const y = cy + radius * Math.sin(-midAngle * RADIAN);

  // Only show label if percentage is greater than 2%
  if (percent < 0.02) return null;

  return (
    <text
      x={x}
      y={y}
      className="text-sm font-medium"
      fill="#374151"
      textAnchor={x > cx ? 'start' : 'end'}
      dominantBaseline="central"
    >
      {`${name} (${(percent * 100).toFixed(1)}%)`}
    </text>
  );
};

const AnalyticsPage = () => {
  const [activeTab, setActiveTab] = useState('amino');
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    fetchAnalyticsData();
  }, []);

  const fetchAnalyticsData = async () => {
    try {
      setLoading(true);
      const analyticsData = await getAnalyticsData();
      setData(analyticsData);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch analytics data');
    } finally {
      setLoading(false);
    }
  };

  const exportData = () => {
    if (!data) return;

    const wb = XLSX.utils.book_new();
    let exportData;
    let sheetName;

    switch(activeTab) {
      case 'amino':
        exportData = data.amino_acids;
        sheetName = 'Amino Acid Composition';
        break;
      case 'length':
        exportData = data.sequence_length;
        sheetName = 'Length Distribution';
        break;
      case 'species':
        exportData = data.species;
        sheetName = 'Species Distribution';
        break;
      case 'functionDescription':
        exportData = data.functionDescription;
        sheetName = 'Function Description';
        break;
      default:
        exportData = [];
        sheetName = 'Data';
    }

    const formattedData = exportData.labels.map((label, index) => ({
      label,
      value: exportData.values[index]
    }));

    const ws = XLSX.utils.json_to_sheet(formattedData);
    XLSX.utils.book_append_sheet(wb, ws, sheetName);
    XLSX.writeFile(wb, `peptide-analytics-${sheetName.toLowerCase()}.xlsx`);
  };

  if (loading) {
    return (
      <div className="min-h-screen flex items-center justify-center bg-gray-50">
        <div className="text-center">
          <Loader2 className="w-12 h-12 animate-spin text-blue-500 mx-auto" />
          <p className="mt-4 text-gray-600">Loading analytics data...</p>
        </div>
      </div>
    );
  }

  if (error || !data) {
    return (
      <div className="min-h-screen flex items-center justify-center bg-gray-50">
        <div className="bg-white rounded-lg shadow-lg p-8 max-w-md w-full text-center">
          <AlertCircle className="w-12 h-12 text-red-500 mx-auto mb-4" />
          <h2 className="text-2xl font-bold text-gray-900 mb-2">Error Loading Data</h2>
          <p className="text-gray-600 mb-6">{error}</p>
          <button
            onClick={fetchAnalyticsData}
            className="px-4 py-2 bg-blue-500 text-white rounded-lg hover:bg-blue-600"
          >
            Retry
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50 p-8">
      <div className="max-w-7xl mx-auto">
        <div className="flex justify-between items-center mb-8">
          <h1 className="text-3xl font-bold text-gray-900">Peptide Analytics</h1>
          <button
            onClick={exportData}
            className="flex items-center px-4 py-2 bg-blue-500 text-white rounded-lg hover:bg-blue-600 transition-colors"
          >
            <Download className="w-4 h-4 mr-2" />
            Export Data
          </button>
        </div>

        <div className="flex space-x-4 mb-8">
          <button
            onClick={() => setActiveTab('amino')}
            className={`px-6 py-2 rounded-lg transition-colors ${
              activeTab === 'amino'
                ? 'bg-blue-500 text-white'
                : 'bg-white text-gray-700 hover:bg-blue-50'
            }`}
          >
            Amino Acid Composition
          </button>
          <button
            onClick={() => setActiveTab('length')}
            className={`px-6 py-2 rounded-lg transition-colors ${
              activeTab === 'length'
                ? 'bg-blue-500 text-white'
                : 'bg-white text-gray-700 hover:bg-blue-50'
            }`}
          >
            Length Distribution
          </button>
          <button
            onClick={() => setActiveTab('species')}
            className={`px-6 py-2 rounded-lg transition-colors ${
              activeTab === 'species'
                ? 'bg-blue-500 text-white'
                : 'bg-white text-gray-700 hover:bg-blue-50'
            }`}
          >
            Species Distribution
          </button>
          <button
            onClick={() => setActiveTab('functionDescription')}
            className={`px-6 py-2 rounded-lg transition-colors ${
              activeTab === 'functionDescription'
                ? 'bg-blue-500 text-white'
                : 'bg-white text-gray-700 hover:bg-blue-50'
            }`}
          >
            Function Description
          </button>
        </div>

        <div className="bg-white rounded-xl shadow-lg p-8">
          {activeTab === 'amino' && (
            <div>
              <h2 className="text-2xl font-bold text-gray-900 mb-6">Amino Acid Composition</h2>
              <div className="h-[600px]">
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart
                    data={data.amino_acids.labels.map((label, index) => ({
                      name: label,
                      count: data.amino_acids.values[index]
                    }))}
                    margin={{ top: 20, right: 30, left: 20, bottom: 5 }}
                  >
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis dataKey="name" />
                    <YAxis />
                    <Tooltip />
                    <Legend />
                    <Bar dataKey="count" name="Frequency" fill="#3B82F6" />
                  </BarChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {activeTab === 'length' && (
            <div>
              <h2 className="text-2xl font-bold text-gray-900 mb-6">Sequence Length Distribution</h2>
              <div className="h-[600px]">
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart
                    data={data.sequence_length.labels.map((label, index) => ({
                      length: label,
                      count: data.sequence_length.values[index]
                    }))}
                    margin={{ top: 20, right: 30, left: 20, bottom: 5 }}
                  >
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis dataKey="length" />
                    <YAxis />
                    <Tooltip />
                    <Legend />
                    <Bar dataKey="count" name="Number of Peptides" fill="#10B981" />
                  </BarChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {activeTab === 'species' && (
            <div>
              <h2 className="text-2xl font-bold text-gray-900 mb-6">Species Distribution</h2>
              <div className="h-[600px]">
                <ResponsiveContainer width="100%" height="100%">
                  <PieChart>
                    <Pie
                      data={data.species.labels
                        .map((label, index) => ({
                          name: label,
                          value: data.species.values[index]
                        }))
                        .sort((a, b) => b.value - a.value)
                        .slice(0, 15)} // Show top 15 species
                      cx="50%"
                      cy="50%"
                      label={CustomPieLabel}
                      labelLine={true}
                      outerRadius={250}
                      fill="#8884d8"
                      dataKey="value"
                    >
                      {data.species.labels.map((_, index) => (
                        <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                      ))}
                    </Pie>
                    <Tooltip />
                  </PieChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {activeTab === 'functionDescription' && (
            <div>
              <h2 className="text-2xl font-bold text-gray-900 mb-6">Function Distribution</h2>
              <div className="h-[600px]">
                <ResponsiveContainer width="100%" height="100%">
                  <PieChart>
                    <Pie
                      data={data.functionDescription.labels
                        .map((label, index) => ({
                          name: label,
                          value: data.functionDescription.values[index]
                        }))
                        .sort((a, b) => b.value - a.value)
                        .slice(0, 12)} // Show top 12 functions
                      cx="50%"
                      cy="50%"
                      label={CustomPieLabel}
                      labelLine={true}
                      outerRadius={250}
                      fill="#8884d8"
                      dataKey="value"
                    >
                      {data.functionDescription.labels.map((_, index) => (
                        <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                      ))}
                    </Pie>
                    <Tooltip />
                  </PieChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default AnalyticsPage;