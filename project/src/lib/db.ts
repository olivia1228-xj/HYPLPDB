export interface PeptideStructure {
  smiles: string;
  inchi: string;
  visualization?: {
    '2d_image': string;
    '3d_structure': {
      atoms: Array<{
        id: number;
        element: string;
        coordinates: {
          x: number;
          y: number;
          z: number;
        };
      }>;
      bonds: Array<{
        start: number;
        end: number;
        order: number;
      }>;
    };
  };
}

export interface ExternalReference {
  id: string;
  link: string;
  display?: string;
}

export interface DatabaseReferences {
  uniprot: ExternalReference;
  biopepuwm: ExternalReference;
  fermfoodb: ExternalReference;
  biopepdb: ExternalReference;
  plantpepdb: ExternalReference;
  cyclicpepedia: ExternalReference;
  jpo: ExternalReference;
}

export interface PeptideExternalReferences {
  literature: ExternalReference[];
  databases: DatabaseReferences;
}

export interface PeptideData {
  basicInfo: PeptideBasicInfo;
  physicalProperties: PeptidePhysicalProperties;
  structure: PeptideStructure;
  externalReferences: PeptideExternalReferences;
}

export interface PeptideBasicInfo {
  id: string;
  species: string;
  sequence: string;
  length: number;
  class: string;
  source: string;
  extractionMethod: string;
  activity: {
    functions: string[];
    functionDescription: string;
    ic50: string;
    validationModel: string;
  };
}

export interface PeptidePhysicalProperties {
  molecularWeight: number;
  isoelectricPoint: number;
  chargeAtPH7: number;
  hydrophobicity: number;
  instabilityIndex: number;
  aromaticity: number;
}

