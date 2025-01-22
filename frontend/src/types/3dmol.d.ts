declare global {
  interface Window {
    $3Dmol: {
      createViewer: (element: HTMLElement, config: any) => any;
      download: (pdb: string, format: string, options?: any) => any;
    }
  }
}

export {};