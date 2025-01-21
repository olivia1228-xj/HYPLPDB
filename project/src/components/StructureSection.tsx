import React, { useEffect, useState } from 'react';
import MoleculeViewer from './MoleculeViewer';

const StructureSection = ({ structure }) => {
  return (
    <div className="space-y-6">
      <h2 className="text-xl font-bold mb-6">Structural Information</h2>

      {/* SMILES and InChI Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
        <div className="bg-gray-50 p-4 rounded-lg">
          <p className="text-sm text-gray-500 mb-2">SMILES</p>
          <p className="font-mono text-sm break-all">{structure.smiles}</p>
        </div>
        <div className="bg-gray-50 p-4 rounded-lg">
          <p className="text-sm text-gray-500 mb-2">InChI</p>
          <p className="font-mono text-sm break-all">{structure.inchi}</p>
        </div>
      </div>

      {/* 2D and 3D Structure Layout */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
        {/* 2D Structure */}
        {structure.visualization?.['2d_image'] && (
          <div className="flex flex-col">
            <h3 className="text-lg font-semibold mb-4">2D Structure</h3>
            <div className="bg-white rounded-lg shadow-md p-4">
              <img
                src={structure.visualization['2d_image']}
                alt="2D structure"
                className="max-w-full h-auto rounded-lg"
              />
            </div>
          </div>
        )}

        {/* 3D Structure */}
        {structure.visualization?.['3d_data'] && (
          <div className="flex flex-col">
            <h3 className="text-lg font-semibold mb-4">3D Structure</h3>
            <div className="bg-white rounded-lg shadow-md">
              <MoleculeViewer pdbData={structure.visualization['3d_data']} />
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default StructureSection;