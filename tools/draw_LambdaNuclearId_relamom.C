#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TString.h"
#include "TStyle.h"
#include <iostream>

void draw_LambdaNuclearId_relamom() {
    // スタイルの設定
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    // 入力ファイル
    TString inputFile = "/star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/auau13p5_anaLambdaNuclearId/auau13p5_anaLambdaNuclearId_alltmp.root";
    TFile* fIn = TFile::Open(inputFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
        return;
    }

    // 出力ディレクトリの作成
    TString outDir = "/star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/auau13p5_anaLambdaNuclearId/pdf";
    gSystem->mkdir(outDir, kTRUE);

    const char* species[] = {"d", "t", "3He", "4He"};
    int nSpecies = 4;
    int nCentBins = 9;

    // 各原子核についてループ
    for (int i = 0; i < nSpecies; ++i) {
        TString sp = species[i];
        
        // 3x3のキャンバスを作成 (9つのCentrality用)
        TCanvas* c1 = new TCanvas(Form("c1_%s", sp.Data()), Form("k* for Lambda - %s", sp.Data()), 1200, 900);
        c1->Divide(3, 3);
        
        // True event histograms are in the "nuclearid" subdirectory
        TDirectory* dirNuc = (TDirectory*)fIn->Get("nuclearid");
        
        for (int j = 0; j < nCentBins; ++j) {
            c1->cd(j + 1);
            
            // ヒストグラムの取得
            TString histName = Form("hKstar_%s_CentBin%d", sp.Data(), j);
            TH1* h = 0;
            if (dirNuc) {
                h = (TH1*)dirNuc->Get(histName);
            }
            if (!h) {
                // fallback: try root directory (for older files without subdirectory structure)
                h = (TH1*)fIn->Get(histName);
            }
            
            if (h) {
                // 統計エラーを考慮 (事前にSumw2が呼ばれていない場合のために)
                if (!h->GetSumw2N()) {
                    h->Sumw2();
                }
                
                // 見やすいように描画スタイルを調整
                h->SetMarkerStyle(20);
                h->SetMarkerSize(0.8);
                h->SetLineColor(kBlack);
                h->SetMarkerColor(kBlack);
                
                // E1 オプションでエラーバー付きで描画
                h->Draw("E1");
            } else {
                std::cerr << "Warning: Histogram " << histName << " not found in file." << std::endl;
            }
        }
        
        // PDFとして保存
        TString outFile = Form("%s/kstar_%s.pdf", outDir.Data(), sp.Data());
        c1->SaveAs(outFile);
        delete c1;
    }
    
    std::cout << "Drawing completed. PDF files saved to: " << outDir << std::endl;
    
    fIn->Close();
}
