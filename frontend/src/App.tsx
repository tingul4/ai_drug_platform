import { useState } from 'react';
import AnalysisPage from './components/AnalysisPage';
import ResultsPage from './components/ResultsPage';

function App() {
  const [activePage, setActivePage] = useState<'analysis' | 'results'>('analysis');
  const [jobId, setJobId] = useState<string | null>(null);

  const handleStartAnalysis = (id: string) => {
    setJobId(id);
    setActivePage('results');
  };

  return (
    <div className="min-h-screen bg-gray-50 text-gray-900 font-sans">
      <header className="bg-white border-b border-gray-200 py-4 px-6 sticky top-0 z-10">
        <div className="max-w-7xl mx-auto flex justify-between items-center">
          <h1 className="text-2xl font-bold text-purple-700 m-0">AIM3 Peptide Platform</h1>
          <nav className="flex space-x-4">
            <button
              onClick={() => setActivePage('analysis')}
              className={`px-3 py-2 rounded-md text-sm font-medium transition-colors ${
                activePage === 'analysis'
                  ? 'bg-purple-100 text-purple-700'
                  : 'text-gray-500 hover:text-gray-700 hover:bg-gray-100'
              }`}
            >
              優化分析
            </button>
            <button
              onClick={() => setActivePage('results')}
              className={`px-3 py-2 rounded-md text-sm font-medium transition-colors ${
                activePage === 'results'
                  ? 'bg-purple-100 text-purple-700'
                  : 'text-gray-500 hover:text-gray-700 hover:bg-gray-100'
              }`}
              disabled={!jobId && activePage !== 'results'}
            >
              分析結果
            </button>
          </nav>
        </div>
      </header>

      <main className="max-w-7xl mx-auto py-8 px-6">
        {activePage === 'analysis' ? (
          <AnalysisPage onStartAnalysis={handleStartAnalysis} />
        ) : (
          <ResultsPage jobId={jobId} onBack={() => setActivePage('analysis')} />
        )}
      </main>
    </div>
  );
}

export default App;
