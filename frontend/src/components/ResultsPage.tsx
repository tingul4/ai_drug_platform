import React, { useState, useEffect } from 'react';
import { ArrowLeft, CheckCircle2, Loader2, BarChart3, TrendingUp, Info, Microscope } from 'lucide-react';
import { getJobStatus, type JobStatus } from '../api';

interface ResultsPageProps {
  jobId: string | null;
  onBack: () => void;
}

const formatDuration = (s: number): string => {
  if (s < 1) return `${s.toFixed(2)}s`;
  if (s < 60) return `${s.toFixed(1)}s`;
  const m = Math.floor(s / 60);
  const sec = Math.round(s - m * 60);
  return `${m}m ${sec}s`;
};

const ResultsPage: React.FC<ResultsPageProps> = ({ jobId, onBack }) => {
  const [jobData, setJobData] = useState<JobStatus | null>(null);
  const [status, setStatus] = useState<'loading' | 'pending' | 'processing' | 'completed' | 'error'>('loading');
  const candidateRows = Array.isArray(jobData?.results?.candidates) ? jobData!.results.candidates : [];
  const alanineRows = Array.isArray(jobData?.results?.alanine_scan) ? jobData!.results.alanine_scan : [];
  const validStatuses = new Set(['pending', 'processing', 'completed', 'error']);

  useEffect(() => {
    if (!jobId) {
      setStatus('error');
      return;
    }

    const pollStatus = async () => {
      try {
        const data = await getJobStatus(jobId) as any;
        if (!data || typeof data !== 'object') {
          console.error('Malformed job payload:', data);
          setStatus('error');
          return true;
        }

        const nextStatus = data.status;
        if (!validStatuses.has(nextStatus)) {
          console.error('Unexpected job status:', nextStatus, data);
          setJobData({ ...data, error: '後端回傳了無效的任務狀態。' } as JobStatus);
          setStatus('error');
          return true;
        }

        if (nextStatus === 'completed' && !data.results) {
          console.error('Completed job missing results:', data);
          setJobData({ ...data, error: '任務已完成，但結果資料缺失。' } as JobStatus);
          setStatus('error');
          return true;
        }

        setJobData(data as JobStatus);
        setStatus(nextStatus);

        if (nextStatus === 'completed' || nextStatus === 'error') {
          return true; // Stop polling
        }
      } catch (error) {
        console.error('Polling error:', error);
        setStatus('error');
        return true;
      }
      return false;
    };

    pollStatus();
    const interval = setInterval(async () => {
      const stop = await pollStatus();
      if (stop) clearInterval(interval);
    }, 2000);

    return () => clearInterval(interval);
  }, [jobId]);

  if (status === 'error') {
    return (
      <div className="text-center py-20 max-w-2xl mx-auto">
        <h2 className="text-2xl font-bold text-red-600">分析出錯</h2>
        {jobData?.error && (
          <pre className="mt-4 text-left text-sm bg-red-50 border border-red-200 rounded-lg p-4 whitespace-pre-wrap text-red-800">
            {jobData.error}
          </pre>
        )}
        {!jobId && <p className="mt-4 text-gray-600">找不到任務 ID</p>}
        <button onClick={onBack} className="mt-4 text-purple-600 hover:underline">返回優化分析</button>
      </div>
    );
  }

  return (
    <div className="space-y-8">
      <div className="flex items-center justify-between">
        <button
          onClick={onBack}
          className="flex items-center space-x-2 text-gray-600 hover:text-purple-600 transition-colors"
        >
          <ArrowLeft className="w-5 h-5" />
          <span>返回設定</span>
        </button>
        <div className="text-sm text-gray-500 text-right space-y-1">
          <div>任務 ID: <code className="bg-gray-100 px-2 py-1 rounded">{jobId}</code></div>
          {jobData?.results?.timings && (
            <div className="text-xs">
              {jobData.results.from_cache ? (
                <span className="text-green-600">⚡ Cache hit ({jobData.results.timings.cache_hit_s}s)</span>
              ) : (
                <span>
                  總時間 <span className="font-medium">{formatDuration(jobData.results.timings.total_s)}</span>
                  {' · '}Tool A {formatDuration(jobData.results.timings.tool_a_s)}
                  {' · '}Tool B {formatDuration(jobData.results.timings.tool_b_s)}
                </span>
              )}
            </div>
          )}
        </div>
      </div>

      {(status === 'processing' || status === 'pending' || status === 'loading') && (
        <section className="bg-white p-12 rounded-xl shadow-sm border border-gray-100 text-center space-y-6">
          <Loader2 className="w-16 h-16 text-purple-600 animate-spin mx-auto" />
          <div className="space-y-2">
            <h2 className="text-2xl font-bold">
              {status === 'pending' ? '任務排隊中...' : '正在進行模擬分析...'}
            </h2>
            <p className="text-gray-500">正在執行 Tool A ({jobData?.params?.tool_a || 'AI'}) 結構生成與 Tool B ({jobData?.params?.tool_b || 'Energy'}) 能量計算</p>
          </div>
          <div className="max-w-md mx-auto">
            <div className="w-full bg-gray-200 rounded-full h-2.5">
              <div 
                className="bg-purple-600 h-2.5 rounded-full transition-all duration-300" 
                style={{ width: `${jobData?.progress || 0}%` }}
              ></div>
            </div>
            <p className="text-right text-sm text-gray-500 mt-2">{jobData?.progress || 0}%</p>
          </div>
        </section>
      )}

      {status === 'completed' && jobData?.results && (
        <div className="space-y-8 animate-in fade-in duration-700">
          {/* Reliability Section */}
          <section className="grid grid-cols-1 md:grid-cols-3 gap-6">
            <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-100 md:col-span-2">
              <div className="flex items-center justify-between mb-6">
                <div className="flex items-center space-x-3">
                  <TrendingUp className="text-green-600 w-6 h-6" />
                  <h3 className="text-xl font-semibold m-0">模擬可靠度分析 (ρ)</h3>
                </div>
                <div className={`px-3 py-1 rounded-full text-sm font-bold ${
                  jobData.results.reliability === 'High' ? 'bg-green-100 text-green-700' :
                  jobData.results.reliability === 'Medium' ? 'bg-yellow-100 text-yellow-700' :
                  'bg-red-100 text-red-700'
                }`}>
                  可信度: {jobData.results.reliability}
                </div>
              </div>
              <div className="flex items-end space-x-8">
                <div className="flex-1">
                  <p className="text-sm text-gray-500 mb-1">Spearman Correlation (ρ)</p>
                  <div className="text-5xl font-black text-purple-700">
                    {jobData.results.rho == null || isNaN(jobData.results.rho)
                      ? '—'
                      : `${jobData.results.rho >= 0 ? '+' : ''}${jobData.results.rho.toFixed(3)}`}
                  </div>
                </div>
                <div className="flex-1 text-sm text-gray-600 border-l border-gray-200 pl-8 space-y-2">
                  <div className="flex justify-between">
                    <span>Tool A:</span>
                    <span className="font-medium">{jobData.results.tool_a}</span>
                  </div>
                  <div className="flex justify-between">
                    <span>Tool B:</span>
                    <span className="font-medium">{jobData.results.tool_b}</span>
                  </div>
                  <div className="flex justify-between">
                    <span>基準數據:</span>
                    <span className="font-medium">Stefansson 1998</span>
                  </div>
                </div>
              </div>
            </div>

            <div className="bg-purple-600 p-6 rounded-xl shadow-md text-white flex flex-col justify-between">
              <div>
                <h3 className="text-lg font-bold mb-2">最佳組合推薦</h3>
                <p className="text-purple-100 text-sm">{jobData.results.tool_a} + {jobData.results.tool_b} 在目前的序列類型中展現了最高的對位精度。</p>
              </div>
              <div className="flex items-center space-x-2 text-sm bg-purple-700 p-3 rounded-lg">
                <CheckCircle2 className="w-4 h-4 text-purple-300" />
                <span>{jobData.results.stefansson_match ? '已重現 Stefansson 熱點方向' : '部分匹配'}</span>
              </div>
            </div>
          </section>

          {/* Pareto Candidates */}
          <section className="bg-white rounded-xl shadow-sm border border-gray-100 overflow-hidden">
            <div className="p-6 border-b border-gray-100 flex items-center justify-between">
              <div className="flex items-center space-x-3">
                <BarChart3 className="text-purple-600 w-6 h-6" />
                <h3 className="text-xl font-semibold m-0">Pareto 候選分子排序 (MPO)</h3>
              </div>
              <button className="text-sm text-purple-600 font-medium hover:underline flex items-center space-x-1">
                <Info className="w-4 h-4" />
                <span>指標說明</span>
              </button>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left">
                <thead className="bg-gray-50 text-gray-500 text-xs uppercase tracking-wider font-semibold">
                  <tr>
                    <th className="px-6 py-4">Rank</th>
                    <th className="px-6 py-4">ID</th>
                    <th className="px-6 py-4">突變點</th>
                    <th className="px-6 py-4">ΔΔG (kcal/mol)</th>
                    <th className="px-6 py-4">免疫風險</th>
                    <th className="px-6 py-4">聚集風險</th>
                    <th className="px-6 py-4">J(x) 綜合評分</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-100">
                  {candidateRows.map((row: any, idx: number) => (
                    <tr key={row.id} className="hover:bg-purple-50/30 transition-colors">
                      <td className="px-6 py-4">
                        <span className={`w-6 h-6 rounded-full flex items-center justify-center text-xs font-bold ${idx === 0 ? 'bg-yellow-100 text-yellow-700' : 'bg-gray-100 text-gray-600'}`}>
                          {idx + 1}
                        </span>
                      </td>
                      <td className="px-6 py-4 font-mono text-sm">{row.id}</td>
                      <td className="px-6 py-4">{row.mutation}</td>
                      <td className="px-6 py-4 text-green-600 font-medium">
                        {row.ddg == null || isNaN(row.ddg) ? '—' : row.ddg}
                      </td>
                      <td className="px-6 py-4">{row.immuno}</td>
                      <td className="px-6 py-4">{row.agg}</td>
                      <td className="px-6 py-4 font-bold text-purple-600">{row.j}</td>
                    </tr>
                  ))}
                  {candidateRows.length === 0 && (
                    <tr>
                      <td className="px-6 py-6 text-sm text-gray-500" colSpan={7}>
                        此組合目前沒有可顯示的候選優化結果。
                      </td>
                    </tr>
                  )}
                </tbody>
              </table>
            </div>
            <div className="p-4 bg-gray-50 border-t border-gray-100 text-center">
              <button className="text-purple-600 text-sm font-semibold hover:text-purple-800">查看完整變體庫</button>
            </div>
          </section>

          {alanineRows.length > 0 && (
            <section className="bg-white rounded-xl shadow-sm border border-gray-100 overflow-hidden">
              <div className="p-6 border-b border-gray-100 flex items-center justify-between">
                <div className="flex items-center space-x-3">
                  <Microscope className="text-purple-600 w-6 h-6" />
                  <h3 className="text-xl font-semibold m-0">Alanine Scan 熱點分析</h3>
                </div>
                <span className="text-xs text-gray-500">
                  {alanineRows.some((r: any) => r.predicted_ddg !== undefined && r.mutant_dG !== undefined)
                    ? '真實 Tool A + Tool B sweep'
                    : '暫用 stub 值（背景 sweep 完成後會替換）'}
                </span>
              </div>
              <div className="overflow-x-auto">
                <table className="w-full text-left">
                  <thead className="bg-gray-50 text-gray-500 text-xs uppercase tracking-wider font-semibold">
                    <tr>
                      <th className="px-6 py-3">位置</th>
                      <th className="px-6 py-3">原始</th>
                      <th className="px-6 py-3">突變</th>
                      <th className="px-6 py-3">ΔΔG (kcal/mol)</th>
                      <th className="px-6 py-3">Kd fold</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-gray-100">
                      {alanineRows.map((row: any) => (
                      <tr key={row.position} className="hover:bg-purple-50/30 transition-colors">
                        <td className="px-6 py-3 font-mono">{row.position}</td>
                        <td className="px-6 py-3 font-mono">{row.original}</td>
                        <td className="px-6 py-3 font-mono">{row.mutant}</td>
                        <td className={`px-6 py-3 font-medium ${
                          row.predicted_ddg != null && !isNaN(row.predicted_ddg)
                            ? row.predicted_ddg > 1.5 ? 'text-red-600' : row.predicted_ddg > 0.5 ? 'text-orange-600' : 'text-gray-700'
                            : 'text-gray-400'
                        }`}>
                          {row.error
                            ? <span className="text-red-500 text-xs font-normal" title={row.error}>計算錯誤 ⚠</span>
                            : row.predicted_ddg == null || isNaN(row.predicted_ddg)
                              ? '—'
                              : row.predicted_ddg.toFixed(2)}
                        </td>
                        <td className="px-6 py-3">
                          {row.kd_fold == null || isNaN(row.kd_fold) ? '—' : row.kd_fold.toFixed(2)}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </section>
          )}
        </div>
      )}
    </div>
  );
};

export default ResultsPage;
