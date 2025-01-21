/// <reference types="vite/client" />

interface Window {
  $3Dmol: {
    createViewer: (element: HTMLElement, config: any) => any;
    download: (type: string, format: string) => void;
  };
  $: any;
}

export {};