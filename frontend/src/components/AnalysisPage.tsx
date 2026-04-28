import React, { useState } from 'react';
import { Play, FileText, Settings } from 'lucide-react';
import { createJob } from '../api';

interface AnalysisPageProps {
  onStartAnalysis: (jobId: string) => void;
}

const AnalysisPage: React.FC<AnalysisPageProps> = ({ onStartAnalysis }) => {
  const [sequence, setSequence] = useState('');
  const [toolA, setToolA] = useState('boltz2');
  const [toolB, setToolB] = useState('mmgbsa');
  const [isSubmitting, setIsSubmitting] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!sequence.trim()) return;

    setIsSubmitting(true);
    try {
      const response = await createJob({ sequence, tool_a: toolA, tool_b: toolB });
      onStartAnalysis(response.job_id);
    } catch (error) {
      console.error('Failed to create job:', error);
      alert('無法建立分析任務，請檢查後端連線。');
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <div className="max-w-3xl mx-auto space-y-8">
      <section className="bg-white p-8 rounded-xl shadow-sm border border-gray-100">
        <div className="flex items-center space-x-3 mb-6">
          <FileText className="text-purple-600 w-6 h-6" />
          <h2 className="text-xl font-semibold m-0">序列輸入與參數設定</h2>
        </div>

        <form onSubmit={handleSubmit} className="space-y-6">
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-2">
              快速選擇測試胜肽 (點擊即可自動填入)
            </label>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
              {[
                { name: 'Peptide A', sequence: 'KGMAPALRHLYKELMGPWNK', desc: 'Wild Type (20 aa)' },
                { name: 'Peptide B', sequence: 'KGMAPALRHLYK', desc: 'Short variant (12 aa)' },
                { name: 'Peptide C', sequence: 'RHLYKELMGPWNK', desc: 'Truncated (13 aa)' }
              ].map((peptide) => (
                <div
                  key={peptide.name}
                  onClick={() => setSequence(peptide.sequence)}
                  className="cursor-pointer border border-gray-200 rounded-lg p-3 hover:border-purple-500 hover:bg-purple-50 transition-colors"
                >
                  <div className="font-semibold text-purple-700">{peptide.name}</div>
                  <div className="text-xs text-gray-800 break-all mt-1 font-mono">{peptide.sequence}</div>
                  <div className="text-xs text-gray-500 mt-1">{peptide.desc}</div>
                </div>
              ))}
            </div>
            
            <label className="block text-sm font-medium text-gray-700 mb-2">
              胜肽胺基酸序列 (FASTA 格式)
            </label>
            <textarea
              value={sequence}
              onChange={(e) => setSequence(e.target.value)}
              placeholder="例如: KGMAPALRHLYK"
              className="w-full h-32 p-4 border border-gray-300 rounded-lg focus:ring-2 focus:ring-purple-500 focus:border-transparent font-mono text-sm"
              required
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Tool A (3D 結構模擬)
              </label>
              <select
                value={toolA}
                onChange={(e) => setToolA(e.target.value)}
                className="w-full p-2.5 border border-gray-300 rounded-lg focus:ring-2 focus:ring-purple-500"
              >
                <option value="boltz2">Boltz-2 (精準度高)</option>
                <option value="colabfold">ColabFold (速度快)</option>
                <option value="chai1">Chai-1</option>
              </select>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Tool B (能量評分/交互作用)
              </label>
              <select
                value={toolB}
                onChange={(e) => setToolB(e.target.value)}
                className="w-full p-2.5 border border-gray-300 rounded-lg focus:ring-2 focus:ring-purple-500"
              >
                <option value="mmgbsa">MM/GBSA (AmberTools, igb=5)</option>
                <option value="iptm">iptm (Tool A confidence)</option>
              </select>
            </div>
          </div>

          <div className="pt-4">
            <button
              type="submit"
              disabled={isSubmitting || !sequence.trim()}
              className="w-full bg-purple-600 hover:bg-purple-700 text-white font-bold py-3 px-6 rounded-lg transition-colors flex items-center justify-center space-x-2 disabled:opacity-50"
            >
              {isSubmitting ? (
                <span>提交中...</span>
              ) : (
                <>
                  <Play className="w-5 h-5" />
                  <span>開始優化分析</span>
                </>
              )}
            </button>
          </div>
        </form>
      </section>

      <section className="bg-purple-50 p-6 rounded-xl border border-purple-100">
        <div className="flex items-center space-x-3 mb-4">
          <Settings className="text-purple-600 w-5 h-5" />
          <h3 className="text-lg font-semibold text-purple-800 m-0">說明與指南</h3>
        </div>
        <ul className="list-disc list-inside text-sm text-purple-700 space-y-2">
          <li>系統將自動對輸入序列進行 Alanine Scanning。</li>
          <li>Tool A 與 Tool B 的組合將影響模擬的可靠度指標 (ρ)。</li>
          <li>分析過程包含 3D 結構生成，可能需要 3-10 分鐘，請耐心等候。</li>
        </ul>
      </section>
    </div>
  );
};

export default AnalysisPage;
