document.addEventListener('DOMContentLoaded', () => {
    // --- DOM Elementleri ---
    const analysisForm = document.getElementById('analysis-form');
    const sequenceInput = document.getElementById('sequence-input');
    const fileInput = document.getElementById('file-input');
    const fileNameSpan = document.getElementById('file-name');
    const sampleButtons = document.querySelectorAll('.sample-btn');
    const analyzeButton = document.getElementById('analyze-button');
    const loadingIndicator = document.getElementById('loading-indicator');
    const errorMessageDiv = document.getElementById('error-message');
    const errorTextSpan = document.getElementById('error-text');

    // Ayar Elementleri
    const orfLengthSlider = document.getElementById('orf-length-slider');
    const orfLengthValueSpan = document.getElementById('orf-length-value');
    const cpgAlgorithmSelect = document.getElementById('cpg-algorithm-select');

    // Sonuç Alanları Container
    const resultsContainer = document.getElementById('results-container');
    const initialInfoDiv = document.getElementById('initial-info');

    // Sonuç Bölümleri
    const summaryResultsDiv = document.getElementById('summary-results');
    const summaryContentDiv = document.getElementById('summary-content');
    const distributionChartCanvas = document.getElementById('distribution-chart'); // Chart için canvas

    const sequenceViewerSection = document.getElementById('sequence-viewer-section');
    const sequenceViewerDiv = document.getElementById('sequence-viewer');

    const visualizationSection = document.getElementById('visualization-section');
    const visualizationDiv = document.getElementById('visualization');
    const annotationDetailsContainer = document.getElementById('annotation-details-container');
    const annotationDetailsDiv = document.getElementById('annotation-details');

    const orfResultsDiv = document.getElementById('orf-results');
    const minOrfLengthSpan = document.getElementById('min-orf-length'); // ORF başlığındaki span
    const orfTableBody = document.getElementById('orf-table-body');
    const noOrfFoundP = document.getElementById('no-orf-found');

    const proteinStatsResultsDiv = document.getElementById('protein-stats-results');
    const proteinStatsContentDiv = document.getElementById('protein-stats-content');
    const proteinHintsSection = document.getElementById('protein-hints-section');
    const proteinHintsList = document.getElementById('protein-hints-list');


    const otherResultsDiv = document.getElementById('other-results'); // Detaylı analizler ana container'ı (details elementi)
    const motifResultsDiv = document.getElementById('motif-results');
    const motifListUl = document.getElementById('motif-list');
    const noMotifFoundP = document.getElementById('no-motif-found'); // Motif bulunamadı mesajı

    const restrictionResultsDiv = document.getElementById('restriction-results');
    const restrictionTableContainer = document.getElementById('restriction-table-container');
    const noRestrictionFoundP = document.getElementById('no-restriction-found');

    const cpgResultsDiv = document.getElementById('cpg-results');
    const cpgListUl = document.getElementById('cpg-list');
    const noCpgFoundP = document.getElementById('no-cpg-found');
    const cpgAlgorithmUsedSpan = document.getElementById('cpg-algorithm-used'); // CpG başlığındaki span

    const framesResultsDiv = document.getElementById('frames-results');
    const framesContentDiv = document.getElementById('frames-content');
    // const noFramesFoundP = document.getElementById('no-frames-found'); // Gerekirse eklenebilir


    // --- Örnek Sekanslar ---
    const sampleSequences = {
        // Gerçek bir insan genomu parçası, TATA kutusu ve ORF başlangıcı içerir
        dna_example: `>NM_078805.4 Drosophila melanogaster DNA replication-related element factor (Dref), transcript variant A, mRNA
ATCAGTGTTAGCGAGAATACTCAACAAATCGCATTTTTTACGACAGTCAGACGTATTGAAATTAAAAAGC
GGTGACTCTCATATTTGAGCAGTCTGTAAAAATAATCTCAAATAAATGAGAGTGCATTTGAAGCATTTAG
CAACGTTTTATTGCCTTATCGATTGACGGCGTGTCGCGTGGCAACTCTGAACTGTACTTCTGATACAGAT
AAAAAACTATTGACTATTGAACCGTCTATCTATAATCGAGAAATTGCGCGCTAATACACGCATTGGGCAC
AGCAATCGTTGTCGAATCTAATTCAAAACAAACAAGAAGATCCCAATCCCCCACAGAAGACAAGATGAGC
GAAGGGGTACCAGCGTCGCCCGTGGCCATTGGCGAACTCAAGTACGAGGATGTGTCGCAGCTAAATTTTT
CCAAACTGTACTCGCCCAAGATGAAAAGCGTATACTGGCGCTACTTTGGATTCCCCTCGAACGACAATAA
TGAGGTGATTACCAAACAGAATGTGGTCTGCATTAAGTGTCACAAGGTGCTGACCAACCATGGCAACACC
ACCAATTTGCGTGCGCATCTCCAGCACCGACACAAGGATCTGTTCAAGGAGCTGTGCCAGGAGCACGACA
TCCATGTGCCGCCGCGCAAGACGCCGCGTAATGTATCCCATCCACCGCTGTCCAAGCGGAATGTTTCCTC
GCGGCGGGTCAAACTGGAGTTCATTAATAACCGGAACCACGATAACGCTTCCGATGATGAACTGGACGAA
GCAGCCGTGGCAACTGCGGCCATGCAGGCGGAGGAGGACGCCTCCTCGCAAACCATGCTGTACGAGACCA
TGGTGCCCTTCACCTACGACGAGGCCGAAAATCTGGTAGAGGAGGAGGAGCGCCTAGTTATGGAACCAAA
GTACGGACGCAAACGCAAAGTGGCCACTCCATCGAGTGCCCTGATGCACGGGCGTGTCATCAAGCACGAG
GAAGGCGGCTATGCCGCCGTGGCGAACATCGCCAATCTGGCCGAAGCCCTTACGGACATTGTGATTAAGG
ATTTGCGCAATGTGGATTCTCTGTACGATGCCGGTTTTGGTGAATTCCTGCGTCAAGTCCTAGGCAATTC
GGCAGCCATGCCCGAGCCACACAAGATCGACTCTTTGATCAATGAGATGCACGCCTCCAAGTTTCTGGAG
ATTGGCGAAATCACGCGCGATTTCACCTCCGAGAAGCCCTTCTCTCTCGCCTTTGAGATGTGGGTGAATG
TTGAACAGCGGCGCTTCCTTAGCATCTTCCATCACTTTCTGGACGAGGAGACGCACTCAGTTCGCGGTAT
GCTGTATGCCACCGTGGAGTACAACGACTACATCGTCTTTGACGATCTGCTGACCGACTTCTATTTGGCC
AATTGCACACTGGCCATCATAAACTACGACGAGGAAGAGGACCTTTTGCACACCTATCTGCGAGAGAAAA
ATATACCCATTTCGCTGTGCTACGTTTCGGTAATTGACAAATGTTTGCGTCGTGTGTTTGAGATCGAAGA
GGTGGCCACTTTATTGGAGCAGGTGAAGGACCTAATGCAGCGTCATTCAACGGAGATTGCCTCCAAGGTG
TCCGAGGTGCCCATGCCCACGTACAACGAGCACTTCCCCTGGACACTGTACGAAACATTGAAGTTCTTTG
CCGAGTCCATATCTTGGTCCGAGGACATGGATCACCTGGTTATATCAGCCAAGACGGTGACGGAGGCCCT
GAGTGCTCTGGTGATTGCTTTGGACACGCTGCGTGGCGAAGACATTCCCTTGTGTAGCATGCTTTCGCCA
ATCACTTCAAAGATTCTTATCAAGAAGTTGGGCATCGCCGAGCAGGACGATCCCCTAATGATGAATTTGA
AGCGGACCATTTCCAGTGTCCTCCAGGCGCACGTTATATCCAATGACAATCTGACCGCTGCTGCATTGTT
GGACCCACGCTTCCACCGCCTGACCACCATTGACAACCTAGAGCGAACCGTTCGTATGCTGACCCATAAG
TATAACATTAACTTTGGTGGAGTGGGCGAAGGGGAATCCAACGAAGTGGCAGCTACCTCCAGCGTGGTGG
CCATTAAGTCAGAACCAAGAGTAGTGGACGGCAGTGCACCAAAGAAGTTGGGTCTGAAACTGCTGTTTGA
CAGTAACGAGATACCAAATCCTCCGAAGCGAGATGCGGACAGCAGCGTGGAGTCCGATCTTAAGCGATAT
CGCAACGAGGTGGTCGTTCAGCTGGATGAGTCGCCCATCGAGTGGTGGCTCAAGATGGGACACATTTATG
GAACGCTCCGCGATTTGGCCAGCCTGTACCATAGTGTGCCCGGCGTGGTGACGCTCAGCTTCAAGAAGGC
GCTGAGAGACCAAATATACGACTTCAACAAGCGATTCATGCTCACCGGTAGTCACATCGACGCCATCCTC
TTTCTTCATCATCACAACAATTAGTCTGCACTGGCGACTACTGTTTCTTATAATGCAAATACTTTTTAAT
TACTTTTAACGTAATCCGATCTTGTGTATATGAAAATGACTACCAGGGGCGAGGCAAATAGTAATAGCGA
TTCCGCAGTACTTTTAACCAAAAGACAGCCTTGGGACGTACACGGAACATCTTAAGTTATCGACGTATAA
ATGTTAAGGGAACATTTTGACATCATATATTCAAATCAAACAGATAATCTAAGAATATTTGTTTAATAAA
TCCATAACGAAAAAATATCACGCTTATCAAGCTCACTTGTAAATCCCCTTATATCCATTTCAAGAGAACT
AATTAAATACTTAACTTAGTTGTACATAATTGATATTTTAATCTTATAAATAAATAAAATGCATCTAACT
TT`,

        // İnsan mRNA dizisinden alınmış, UTR ve ORF içeren örnek (örnek olarak GAPDH)
        rna_example: `>LT883170.1 Nicotiana benthamiana mRNA for RNA (NbPPA1.1 gene), strain PRJNA188486
ATGGGAGCGCAAATTTTATCAGATCTCGGAACGGAGATTTTGATTCCGGTATGTGCCGTCGTCGGCATCG
CTTTCTCGCTTTTTCAATGGTTTCTCGTCTCTAAAGTGACGCTCAGTGCTGAAAAGTCCTCCGGCGCCGC
CGACGATAAGAATGGTTACGCTGCTGAATCGCTCATTGAAGAAGAAGAAGGCATTAATGACCATAGCGTT
GTTCAGAAATGCGCCGAAATTCAGAACGCCATTTCCGAAGGTGCAACATCATTCCTTTTTACTGAGTACC
AGTATGTTGGTGTTTTCATGGTTGCTTTTGCAATTTTAATCTTCCTTTTCCTCGGCTCTGTCGAGGGTTT
CAGCACAAAGAACCAATCCTGTACATATGACACCACCAAAACTTGCAAGCCCGCTCTTGCCACTGCTGTC
TTCAGCACTGTATCTTTCTTGCTTGGTGCTGTTACATCTGTGGTGTCCGGGTTCCTTGGAATGAAAATTG
CCACATATGCAAATGCCAGAACTACCTTGGAGGCTAGAAAAGGTGTTGGGAAGGCTTTTATTGTTGCATT
CAGATCTGGTGCTGTAATGGGTTTCCTTCTTGCTGCAAACGGTCTTCTGGTTTTGTACATAACCATCCTT
CTATTTAAGTTGTACTATGGTGATGACTGGGAAGGTCTTTTTGAGGCTATAACAGGTTATGGGCTTGGTG
GATCTTCAATGGCCCTTTTTGGTAGAGTTGCTGGAGGTATTTATACCAAAGCAGCAGATGTTGGAGCTGA
TCTTGTGGGCAAGGTGGAAAGGAACATCCCTGAAGATGATCCTAGAAACCCAGCGGTTATTGCTGACAAT
GTCGGTGACAATGTCGGAGATATTGCTGGTATGGGATCAGATCTGTTTGGGTCTTATGCAGAGTCATCCT
GTGCAGCTCTTGTTGTTGCTTCAATCTCCTCCTTCGGTGTCAACCATGAGTTCACTGCTATGCTATATCC
CCTTCTCGTGAGCTCTGTTGGTATCCTCGTTTGTTTGCTTACCACCTTATTCGCAACTGATTTCTTTGAA
GTCAAGGCTGTTAAGGAAATCGAGCCAGCATTGAAGCAGCAACTCGTTATCTCAACTGCTCTCATGACAG
TTGGAATTGCTGCTGTTACTTGGATTGCCCTTCCATCGACCTTCACAATATTCAATTTTGGGGCTCAGAA
AGAAGTAAAAAGCTGGCAGTTGTTCTTGTGTGTCGGAGTTGGTTTGTGGGCTGGACTTATTATTGGGTTC
GTCACCGAGTACTATACCAGCAATGCTTACAGCCCTGTGCAAGATGTTGCTGATTCATGTCGGACCGGTG
CTGCAACAAATGTTATTTTTGGCCTTGCCTTGGGTTATAAATCAGTGATCATCCCCATATTTGCCATAGC
AGTCAGTATATTTGTTAGCTTCAGCTTTGCTGCAATGTATGGCATTGCAGTTGCTGCCCTAGGAATGCTG
AGCACTATAGCCACTGGTTTGGCAATTGATGCATATGGTCCCATCAGTGATAATGCAGGAGGCATTGCTG
AGATGGCTGGAATGAGCCAGAGAATCAGAGAGAGAACTGATGCACTTGATGCTGCAGGAAACACCACTGC
AGCTATTGGAAAGGGATTTGCCATTGGTTCTGCTGCTCTCGTGTCTCTGGCTCTCTTTGGGGCATTTGTC
AGCCGAGCAGCAATTTCCACTGTAGATGTTTTGACGCCTAAAGTATTTATTGGTTTGCTAGTCGGTGCCA
TGCTTCCTTATTGGTTCTCCGCCATGACAATGAAGAGTGTTGGAAGCGCTGCTCTTAAGATGGTTGAGGA
AGTGCGCAGGCAATTCAACACCATTCCTGGTCTCATGGAAGGAACTGCCAAGCCTGACTACGCCACCTGC
GTCAAGATCTCTACAGATGCATCGATCAAGGAGATGATTGCACCTGGTGCTCTTGTCATGCTCACTCCAT
TGATTGTTGGAATCTTGTTTGGCGTCGAAACACTTTCTGGTGTGCTTGCAGGATCTCTCGTCTCTGGTGT
ACAGATTGCCATCTCTGCATCTAATACAGGTGGTGCCTGGGATAATGCTAAGAAGTACATTGAGGCCGGA
GTCTCAGAACATGCAAGAACCCTTGGCCCCAAGGGATCTGATGCACACAAGGCTGCTGTTATTGGTGATA
CTGTTGGTGACCCTCTCAAGGACACATCTGGACCATCATTGAACATTCTTATCAAGCTGATGGCTGTCGA
ATCCCTCGTGTTTGCTCCCTTCTTTGCAACTCATGGTGGTCTTCTCTTCAAGCTATTTTAA
`,

        // Gerçek bir protein dizisi: İnsan İnsülin (NCBI: NP_000198.1)
        protein_example: `>KAI4357646.1 hypothetical protein L6164_001582 [Bauhinia variegata]
MAESVVSFLLNELSKVLKEEINLERGVKREIELIKAELGSINAFLKVADVEEEKDPGLAEWVKQVRDVAF
DIEDILDTFLLHLVARYPRFKCLNKIASSYKKWKINRQTGSELQSIKTKVKEIYERHQRYPSLLEPGLST
RNTKMEEVNQFRRNALLQDKSQLVGTEEPKQELIDWLMQGDSHLKVLAISGMGGLGKTTLVQQIYSSDKV
NARFESKALVTVSQTLKQEELFRELIRQLSPPELHHGIDAMDNHRWKQTTKELLQDKKYLIVFDDVWATQ
DWDQLIVAMPTNSHGSRVIITTRISEVARYSTNVNQGMSYPMKVLSESDSWTLFCKKTFQAESCPPHLKD
ICEQIIKKCGGVPLAVSTISGLLATKGVESFGVWKMILDSLSVELQENSRLQMMQRIFSLSFHDLPYNLK
SCFLYLSIFPEDYLIERTRVVRLWIAEGFVKVKDGRTLEEVAESYFYELLNRSLIQAATTTYDGRVKSCR
IHDFLREAILSKSKDENFVQVLSRENATRLPEGTRRLSIFHLSQQVAMYQFKRSLLSKLRSLFFFGADSS
LTMSELLSGGFKLLKVLDLQGAPLDKCPKQVSTLFLLRYLSFRNTKVSMLSKSIGNLHSLETLDLRQTNV
TELPSEIYKLRKLRNLLICQFHSSGTKCALKVVKSGGCFVGYFIWSVCAFEASKEIGGLVSLQKLCYVDA
SERDDLIIELGKLVQLRRLCLINVKGQHGSVLCSSISKMQCLSSLSLHWKGPIDLQSMSSPPPLLKRLHL
TGRLVELPSWLSSLQYLEVLHLVNSELKEDPLKSLPPGGCLIDLLLYDAFIGEVMFIQPGGCFTRLRFLG
IFKFGSVSKLLIFPNAIPSLECVYIDGEDTAYDGRKGYIDDGEMSDWMPEKGWYQYAEKEFAHIVKPRT`,
    };


    // --- Global Değişkenler ---
    let currentAnnotations = []; // Görselleştirme ve detaylar için
    let currentSequenceType = null; // Mevcut sekans türünü sakla
    let hasOtherResults = false; // Detaylı analiz sonuçlarının varlığını kontrol etmek için

    // --- Yardımcı Fonksiyonlar ---
    function showLoading() {
        analyzeButton.disabled = true;
        loadingIndicator.classList.remove('hidden');
        errorMessageDiv.classList.add('hidden'); // Önceki hatayı gizle
        resultsContainer.classList.add('hidden-initially'); // Sonuçları hemen gizle (opacity 0)
        resultsContainer.classList.remove('visible-now');   // Görünürlük sınıfını kaldır
        initialInfoDiv.classList.add('hidden'); // Başlangıç mesajını gizle
        clearResults(); // Sonuçları temizle (içerikleri boşalt)
    }

    function hideLoading() {
        analyzeButton.disabled = false;
        loadingIndicator.classList.add('hidden');
    }

    function showError(message) {
        if (errorTextSpan) errorTextSpan.textContent = message;
        if (errorMessageDiv) errorMessageDiv.classList.remove('hidden');
        if (resultsContainer) resultsContainer.classList.add('hidden-initially'); // Hata varsa sonuçları gizli tut
        if (initialInfoDiv) initialInfoDiv.classList.remove('hidden'); // Başlangıç mesajını tekrar göster
    }

    function clearResults() {
        // Ana sonuç container'ını gizle (ama DOM'dan kaldırma)
        if (resultsContainer) {
            resultsContainer.classList.add('hidden-initially');
            resultsContainer.classList.remove('visible-now');
        }

        // Tek tek tüm sonuç BÖLÜMLERİNİ gizle
        const sectionsToHide = [
            summaryResultsDiv, sequenceViewerSection, visualizationSection,
            orfResultsDiv, proteinStatsResultsDiv, otherResultsDiv, motifResultsDiv,
            restrictionResultsDiv, cpgResultsDiv, framesResultsDiv, annotationDetailsContainer
        ];
        sectionsToHide.forEach(section => section?.classList.add('hidden'));

        // İçeriklerini temizle
        const contentsToClear = [
            summaryContentDiv, sequenceViewerDiv, visualizationDiv,
            orfTableBody, proteinStatsContentDiv, motifListUl,
            restrictionTableContainer, cpgListUl, framesContentDiv,
            annotationDetailsDiv, proteinHintsList
        ];
        contentsToClear.forEach(content => { if (content) content.innerHTML = ''; });

        // "Bulunamadı" mesajlarını gizle
        const noDataMessages = [
            noOrfFoundP, noMotifFoundP, noRestrictionFoundP, noCpgFoundP
        ];
        noDataMessages.forEach(msg => msg?.classList.add('hidden'));

        // Grafik alanını temizle
        if (window.distributionChart) {
            window.distributionChart.destroy();
            window.distributionChart = null;
        }
        if (distributionChartCanvas && distributionChartCanvas.parentElement) {
            // Canvas'ı yeniden oluşturmak yerine, sadece container'ı temizleyebiliriz
            // Ancak Chart.js destroy sonrası canvas'ı yeniden kullanabilir. Şimdilik sadece chart objesini null yapalım.
            // distributionChartCanvas.parentElement.innerHTML = '<canvas id="distribution-chart"></canvas>'; // Canvas'ı yeniden ekle
        }


        // Görselleştirme alanını başlangıç durumuna getir
        if (visualizationDiv) visualizationDiv.innerHTML = '<span class="text-gray-400 text-sm italic">Analiz bekleniyor...</span>';
        if (annotationDetailsContainer) annotationDetailsContainer.classList.add('hidden');

        // Anotasyonları temizle
        currentAnnotations = [];

        // Hata mesajını gizle
        if (errorMessageDiv) errorMessageDiv.classList.add('hidden');
    }

    // DNA/RNA/Protein baz/AA renklendirme (Sequence Viewer & Visualization için)
    function getBaseClass(base, type) {
        base = base.toUpperCase();
        if (type === 'DNA') {
            switch (base) {
                case 'A': return 'base-A'; case 'T': return 'base-T';
                case 'G': return 'base-G'; case 'C': return 'base-C';
                default: return 'base-Unknown'; // N vb. için
            }
        } else if (type === 'RNA') {
            switch (base) {
                case 'A': return 'base-A'; case 'U': return 'base-U';
                case 'G': return 'base-G'; case 'C': return 'base-C';
                default: return 'base-Unknown';
            }
        } else if (type === 'Protein') {
            // Amino asit gruplarına göre sınıf döndür
            if (['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'G'].includes(base)) return 'amino-acid-hydrophobic'; // Geniş hidrofobik + Glisin/Prolin
            if (['S', 'T', 'C', 'Y', 'N', 'Q'].includes(base)) return 'amino-acid-polar';
            if (['K', 'R', 'H'].includes(base)) return 'amino-acid-positive';
            if (['D', 'E'].includes(base)) return 'amino-acid-negative';
            //if (['C', 'G', 'P'].includes(base)) return 'amino-acid-special'; // Özel grup ayrımı yerine hidrofobik/polar içine dahil edilebilir
            return 'amino-acid-unknown'; // X, *, - vb.
        } else {
            return 'base-Unknown';
        }
    }

    // Görselleştirme Alanını Oluşturma
    function renderVisualization(sequence, annotations, type) {
        if (!visualizationDiv || !sequence) {
            console.error("Visualization div or sequence not found!");
            if (visualizationDiv) visualizationDiv.innerHTML = '<span class="text-red-500">Görselleştirme oluşturulamadı.</span>';
            return;
        }
        visualizationDiv.innerHTML = ''; // Temizle
        if (annotationDetailsContainer) annotationDetailsContainer.classList.add('hidden'); // Detayları gizle
        currentAnnotations = annotations || []; // Global anotasyonları güncelle

        const fragment = document.createDocumentFragment();
        const sequenceLength = sequence.length;

        // Özelliklerin pozisyonlara göre haritasını çıkar (performans için)
        // Key: position, Value: array of annotation objects starting/covering this position
        const featuresMap = new Map();
        currentAnnotations.forEach((anno, index) => {
            if (typeof anno.position !== 'number' || anno.position < 0) return; // Geçersiz pozisyonu atla

            const startPos = anno.position;
            // Bitiş pozisyonunu belirle (dahil - inclusive)
            let endPos = anno.end_position !== undefined ? anno.end_position : startPos;
            // Eğer end_position başlangıçtan küçükse (örn. ters strand ORF), düzelt
            // Ters strand ORF'larda position > end_position olabilir. Bu durumda sadece position'ı işaretlemek yerine
            // aralığı doğru şekilde ele almak gerekir. Backend'den gelen veride end_position orf'un son bazının indexi (dahil) olmalı.
            // Eğer endPos < startPos geliyorsa, bu muhtemelen ters strand ve aralığı doğru hesaplamamız gerekiyor.
            // Basit bir düzeltme: Eğer startPos ve endPos arasındaki fark çok büyükse veya negatifse, bu aralığı
            // olduğu gibi kullanmak yerine sadece startPos'u işaretlemek daha güvenli olabilir
            // Şimdilik backend'in doğru end_position (dahil) döndürdüğünü varsayalım. Eğer start > end ise swap yapalım.
            if (endPos < startPos) {
                [startPos, endPos] = [endPos, startPos]; // Swap
            }


            // Ensure endPos is within bounds
            endPos = Math.min(endPos, sequenceLength - 1);


            // Özelliğin kapladığı her pozisyona işaret koy
            for (let i = startPos; i <= endPos; i++) {
                if (i < 0 || i >= sequenceLength) continue; // Sınır kontrolü
                if (!featuresMap.has(i)) {
                    featuresMap.set(i, []);
                }
                // Orijinal index'i ve tipi sakla
                featuresMap.get(i).push({ ...anno, originalIndex: index });
            }
        });

        // Satır uzunluğu (her satırda kaç karakter olacak)
        const lineLength = 60;
        let positionCounter = 0;

        // Her satırı ayrı ayrı oluştur
        for (let lineStart = 0; lineStart < sequenceLength; lineStart += lineLength) {
            // Satır başına pozisyon numarası ekle
            const positionSpan = document.createElement('span');
            positionSpan.className = 'position-number';
            positionSpan.textContent = (positionCounter + 1).toString();
            fragment.appendChild(positionSpan);
            fragment.appendChild(document.createTextNode(' '));

            // Satırdaki her baz için
            const lineEnd = Math.min(lineStart + lineLength, sequenceLength);
            for (let i = lineStart; i < lineEnd; i++) {
                const base = sequence[i];
                const span = document.createElement('span');
                span.textContent = base;
                // Temel renk sınıfını al
                span.className = `base ${getBaseClass(base, type)}`; // Renk sınıfı
                span.dataset.position = i; // 0-tabanlı pozisyon

                // Bu pozisyonda özellik var mı?
                if (featuresMap.has(i)) {
                    const featuresAtPos = featuresMap.get(i);
                    // Özellikleri önceliğe göre sırala (örn: ORF > codon > motif) - isteğe bağlı
                    featuresAtPos.sort((a, b) => {
                        const priority = { 'ORF': 1, 'cpg_island': 2, 'start_codon': 3, 'stop_codon': 3, 'restriction_site': 4, 'motif': 5 };
                        const typeA = a.type;
                        const typeB = b.type;
                        return (priority[typeA] || 99) - (priority[typeB] || 99);
                    });

                    const mainFeature = featuresAtPos[0]; // En öncelikli özelliği al

                    // Stil sınıfını ekle
                    span.classList.add('feature'); // Genel özellik işareti
                    let featureClass = '';
                    if (mainFeature.type === 'start_codon') featureClass = 'start-codon';
                    else if (mainFeature.type === 'stop_codon') featureClass = 'stop-codon';
                    else if (mainFeature.type === 'ORF') featureClass = 'orf';
                    else if (mainFeature.type === 'restriction_site') featureClass = 'restriction-site';
                    else if (mainFeature.type === 'cpg_island') featureClass = 'cpg-island';
                    // Motifler için genel bir sınıf veya alt tiplere göre sınıflar
                    else if (mainFeature.type?.includes('motif') || mainFeature.type?.includes('box') || mainFeature.type?.includes('site') || mainFeature.type?.includes('element') || mainFeature.type === 'CpG') featureClass = 'motif'; // Genel motif sınıfı

                    if (featureClass) span.classList.add(featureClass);

                    // Tooltip için açıklama ekle
                    if (mainFeature.desc) {
                        span.title = mainFeature.desc;
                    } else {
                        // Eğer açıklama yoksa, basit bir tooltip oluştur
                        let tooltipText = `${mainFeature.type} (Poz: ${mainFeature.position + 1}`;
                        if (mainFeature.end_position !== undefined && mainFeature.end_position !== mainFeature.position) {
                            tooltipText += `-${mainFeature.end_position + 1}`;
                        }
                        if (mainFeature.enzyme) tooltipText += `, Enzim: ${mainFeature.enzyme}`;
                        if (mainFeature.frame) tooltipText += `, Çerçeve: ${mainFeature.frame}`;
                        tooltipText += ')';
                        if (mainFeature.sequence) tooltipText += `\nDizi: ${mainFeature.sequence}`;

                        span.title = tooltipText;
                    }

                    // Tıklama/Hover için index
                    span.dataset.annotationIndex = mainFeature.originalIndex;
                }
                fragment.appendChild(span);
            }

            // Satır sonuna pozisyon numarası ekle
            fragment.appendChild(document.createTextNode(' '));
            const endPositionSpan = document.createElement('span');
            endPositionSpan.className = 'position-number';
            endPositionSpan.textContent = (positionCounter + (lineEnd - lineStart)).toString();
            fragment.appendChild(endPositionSpan);

            // Satır sonu
            fragment.appendChild(document.createElement('br'));
            positionCounter += (lineEnd - lineStart);
        }

        visualizationDiv.appendChild(fragment);
        if (visualizationSection) visualizationSection.classList.remove('hidden'); // Bölümü göster

        // Görselleştirme alanına tıklama olayı ekle
        visualizationDiv.addEventListener('click', function (e) {
            const target = e.target.closest('span.feature'); // Sadece feature class'ına sahip spanları hedefle
            if (target && target.dataset.annotationIndex) {
                const index = parseInt(target.dataset.annotationIndex);
                if (!isNaN(index) && index >= 0 && index < currentAnnotations.length) {
                    showAnnotationDetails(currentAnnotations[index]);
                }
            }
        });

        // Mouseleave olayı ekleyerek alan dışına çıkınca detayları gizle
        visualizationDiv.addEventListener('mouseleave', function (e) {
            // Eğer mouse visualizationDiv dışına çıkıyorsa detayları gizle
            // event.relatedTarget null ise pencere dışına çıkılmıştır
            if (!e.relatedTarget || !visualizationDiv.contains(e.relatedTarget)) {
                hideAnnotationDetails();
            }
        });

        // Mouseover olayı ekleyerek feature üzerine gelince detayları göster
        visualizationDiv.addEventListener('mouseover', (event) => {
            const target = event.target.closest('span.feature');
            if (target && target.dataset.annotationIndex) {
                const index = parseInt(target.dataset.annotationIndex);
                if (!isNaN(index) && currentAnnotations && currentAnnotations[index]) {
                    showAnnotationDetails(currentAnnotations[index]);
                }
            } else {
                // Eğer fare bir feature üzerinde değilse (ama visualizationDiv içindeyse), detayları gizleme
                // (Tooltip gibi davranmasını istemiyorsak)
                // hideAnnotationDetails(); // Tooltip gibi davranması isteniyorsa burayı açın
            }
        });

    }

    // Sekans Görüntüleyici (Sadece Renkli Sekans)
    function renderSequenceViewer(sequence, type) {
        if (!sequenceViewerDiv || !sequence) {
            console.error("Sequence viewer div or sequence not found!");
            if (sequenceViewerDiv) sequenceViewerDiv.innerHTML = '<span class="text-red-500">Sekans görüntülenemedi.</span>';
            return;
        }
        sequenceViewerDiv.innerHTML = ''; // Temizle

        const fragment = document.createDocumentFragment();
        for (let i = 0; i < sequence.length; i++) {
            const base = sequence[i];
            const nucleotideSpan = document.createElement('span');
            nucleotideSpan.textContent = base;
            // Temel renk/stil sınıfını al
            nucleotideSpan.className = `nucleotide ${getBaseClass(base, type)}`;
            fragment.appendChild(nucleotideSpan);
        }
        sequenceViewerDiv.appendChild(fragment);
        if (sequenceViewerSection) sequenceViewerSection.classList.remove('hidden');
    }

    // Protein İstatistiklerini ve Formatlı Diziyi Görüntüleme
    function displayProteinStats(statsData) {
        if (!proteinStatsResultsDiv || !proteinStatsContentDiv) return;

        if (statsData && !statsData.error) {
            proteinStatsContentDiv.innerHTML = `
                <div class="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-2 text-sm">
                    <div><span class="font-medium text-gray-600">Uzunluk:</span> ${statsData.length || 'N/A'} aa</div>
                    <div><span class="font-medium text-gray-600">Molekül Ağırlığı:</span> ${statsData.molecular_weight !== undefined ? statsData.molecular_weight.toLocaleString('tr-TR') + ' Da' : 'N/A'}</div>
                    <div><span class="font-medium text-gray-600">İzoelektrik Nokta (pI):</span> ${statsData.isoelectric_point !== undefined ? statsData.isoelectric_point : 'N/A'}</div>
                    <div><span class="font-medium text-gray-600">Aromatiklik:</span> ${statsData.aromaticity !== undefined ? statsData.aromaticity : 'N/A'}</div>
                    <div><span class="font-medium text-gray-600">Kararlılık İndeksi:</span> ${statsData.instability_index !== undefined ? statsData.instability_index + ' (' + statsData.stability_class + ')' : 'N/A'}</div>
                    <div><span class="font-medium text-gray-600">Hidrofobisite (GRAVY):</span> ${statsData.gravy !== undefined ? statsData.gravy + ' (' + statsData.hydropathy_class + ')' : 'N/A'}</div>
                </div>
                <div class="mt-4 text-sm">
                    <span class="font-medium text-gray-600">İkincil Yapı Tahmini (%):</span>
                    <span class="ml-2">Helix: ${statsData.secondary_structure_fraction?.Helix ?? 'N/A'}%</span>,
                    <span class="ml-2">Turn: ${statsData.secondary_structure_fraction?.Turn ?? 'N/A'}%</span>,
                    <span class="ml-2">Sheet: ${statsData.secondary_structure_fraction?.Sheet ?? 'N/A'}%</span>
                </div>
                 <div class="mt-4 text-sm">
                     <details>
                         <summary class="text-sm font-medium text-gray-600 cursor-pointer hover:text-indigo-600">Amino Asit Grupları (%)</summary>
                         <div class="text-xs mt-1 pl-4">
                            ${Object.entries(statsData.amino_acid_groups || {}).map(([group, perc]) => `<div><span class="capitalize">${group.replace('_', ' ')}:</span> ${perc}%</div>`).join('')}
                         </div>
                     </details>
                 </div>
                <div class="mt-4 text-sm">
                    <details>
                         <summary class="text-sm font-medium text-gray-600 cursor-pointer hover:text-indigo-600">Formatlanmış Protein Dizisi</summary>
                         <pre class="mt-2 text-[10px] md:text-xs font-mono bg-gray-50 p-3 border border-gray-200 rounded max-h-80 overflow-auto protein-sequence-container">${statsData.formatted_sequence || 'N/A'}</pre>
                    </details>
                </div>
                 <div class="mt-4 text-sm">
                     <details>
                         <summary class="text-sm font-medium text-gray-600 cursor-pointer hover:text-indigo-600">Amino Asit Yüzdeleri</summary>
                         <div class="text-xs mt-1 pl-4 grid grid-cols-2 sm:grid-cols-3 md:grid-cols-4 gap-x-3 gap-y-1">
                             ${Object.entries(statsData.amino_acid_percent_named || {}).sort((a, b) => a[0].localeCompare(b[0])).map(([name, perc]) => `<div>${name}: ${perc}%</div>`).join('')}
                         </div>
                     </details>
                 </div>
            `;

            // Eğitimsel Fonksiyon İpuçları (Varsa)
            const hints = statsData.potential_function_hints;
            if (proteinHintsSection && proteinHintsList && hints && hints.length > 0) {
                proteinHintsList.innerHTML = hints.map(hint => `<li>${hint}</li>`).join('');
                proteinHintsSection.classList.remove('hidden');
            } else if (proteinHintsSection) {
                proteinHintsSection.classList.add('hidden');
            }

            proteinStatsResultsDiv.classList.remove('hidden');

        } else if (statsData && statsData.error) {
            proteinStatsContentDiv.innerHTML = `<p class="text-red-500 italic text-sm">Protein istatistikleri hesaplanamadı: ${statsData.error}</p>`;
            proteinStatsResultsDiv.classList.remove('hidden'); // Hata mesajını göstermek için
            if (proteinHintsSection) proteinHintsSection.classList.add('hidden'); // İpuçları bölümünü gizle
        } else {
            // Veri yoksa gizle
            proteinStatsResultsDiv.classList.add('hidden');
        }
    }


    // Ana Sonuç Gösterim Fonksiyonu
    function displayResults(data) {
        try {
            // --- Temizlik ve Hazırlık ---
            clearResults(); // Önce temizle
            if (!resultsContainer) {
                console.error("Results container not found!");
                showError("Arayüz hatası: Sonuç container elementi bulunamadı.");
                return;
            }
            resultsContainer.classList.remove('hidden-initially'); // Ana container'ı görünür yapmaya hazırlan (opacity için)
            initialInfoDiv.classList.add('hidden'); // Başlangıç mesajını gizle
            currentSequenceType = data.type; // Sekans türünü sakla

            // --- 1. Özeti Göster ---
            if (summaryResultsDiv && summaryContentDiv) {
                summaryContentDiv.innerHTML = `
                    <div><span class="font-medium text-gray-600">Sekans Türü:</span> ${data.type || 'Bilinmiyor'}</div>
                    <div><span class="font-medium text-gray-600">Uzunluk:</span> ${data.length || 0} ${data.type === 'Protein' ? 'aa' : 'bp'}</div>
                    ${data.gc_content !== null && data.gc_content !== undefined ? `<div><span class="font-medium text-gray-600">GC İçeriği:</span> ${data.gc_content}%</div>` : ''}
                    ${data.type !== 'Protein' ? `<div><span class="font-medium text-gray-600">Başlangıç Kodonları:</span> ${data.start_codons ? data.start_codons.length : 0}</div>` : ''}
                    ${data.type !== 'Protein' ? `<div><span class="font-medium text-gray-600">Bitiş Kodonları:</span> ${data.stop_codons ? data.stop_codons.length : 0}</div>` : ''}
                    ${data.type !== 'Protein' ? `<div><span class="font-medium text-gray-600">Bulunan ORF Sayısı:</span> ${data.orfs ? data.orfs.length : 0} ${data.orfs && data.orfs.length > 0 ? '(En uzunu: ' + data.orfs[0].length_aa + ' aa)' : ''}</div>` : ''}
                    ${data.type !== 'Protein' ? `<div><span class="font-medium text-gray-600">Bulunan Motifler:</span> ${data.motifs ? data.motifs.length : 0}</div>` : ''}
                     ${data.type === 'DNA' ? `
                    <div><span class="font-medium text-gray-600">CpG Adası Sayısı:</span> ${data.cpg_islands ? data.cpg_islands.length : 0}</div>
                    <div><span class="font-medium text-gray-600">Restriksiyon Alanı Sayısı (Enzimler):</span> ${data.restriction_sites && typeof data.restriction_sites === 'object' && !data.restriction_sites.error ? Object.keys(data.restriction_sites).length : 0}</div>
                     ` : ''}
                 `;
                summaryResultsDiv.classList.remove('hidden');

                // Dağılım grafiğini oluştur (Canvas hazırsa)
                // Dağılım verisi sadece DNA/RNA içinse GC bazlı renk, proteinse genel renk kullansın
                if (distributionChartCanvas && data.distribution && Object.keys(data.distribution).length > 0 && typeof window.createDistributionChart === 'function') {
                    window.createDistributionChart(data.distribution, data.type); // Sekans tipini de pass et
                } else if (distributionChartCanvas && distributionChartCanvas.parentElement) {
                    // Veri yoksa veya hata varsa grafik alanına mesaj yaz
                    distributionChartCanvas.parentElement.innerHTML = '<p class="text-center text-xs text-gray-500 italic p-4">Dağılım grafiği için veri yok veya yetersiz.</p>';
                }
            }

            // --- 2. Sekans Görüntüleyiciyi Oluştur ---
            // Visualization bölümü sequence viewer'ın yerini alabilir veya ayrı tutulabilir
            // Şimdilik visualisation'ı daha detaylı olan olarak kullanalım.
            // if (data.cleaned_sequence && sequenceViewerDiv && sequenceViewerSection) {
            //     renderSequenceViewer(data.cleaned_sequence, data.type);
            // }


            // --- 3. İnteraktif Görselleştirmeyi Oluştur ---
            if (data.cleaned_sequence && data.annotations && visualizationDiv && visualizationSection) {
                // Çok uzun sekanslar için görselleştirmeyi limitle? (Örn. ilk 10000 baz/aa)
                const MAX_VIS_LENGTH = 10000; // Örnek limit
                let sequenceToVisualize = data.cleaned_sequence;
                let annotationsToVisualize = data.annotations;
                let visMessage = '';

                if (data.cleaned_sequence.length > MAX_VIS_LENGTH) {
                    sequenceToVisualize = data.cleaned_sequence.substring(0, MAX_VIS_LENGTH);
                    annotationsToVisualize = data.annotations.filter(a => a.position < MAX_VIS_LENGTH);
                    visMessage = `<p class="text-orange-600 text-sm p-4 italic">Sekans çok uzun (${data.cleaned_sequence.length} ${data.type === 'Protein' ? 'aa' : 'bp'}), performans nedeniyle interaktif görselleştirmenin yalnızca ilk ${MAX_VIS_LENGTH} ${data.type === 'Protein' ? 'aa' : 'bazı'} gösteriliyor.</p>`;
                    visualizationDiv.innerHTML = visMessage; // Mesajı visualizationDiv'e ekle
                } else {
                    visualizationDiv.innerHTML = ''; // Mesaj alanı boşsa temizle
                }


                renderVisualization(sequenceToVisualize, annotationsToVisualize, data.type);

            } else if (visualizationDiv) {
                visualizationDiv.innerHTML = '<span class="text-red-500 italic">Görselleştirilecek sekans veya anotasyon verisi alınamadı.</span>';
                if (visualizationSection) visualizationSection.classList.remove('hidden');
            }


            // --- 4. ORF Tablosunu Doldur ---
            if (orfResultsDiv && orfTableBody && noOrfFoundP && minOrfLengthSpan) {
                orfTableBody.innerHTML = ''; // Temizle
                if (data.type === 'DNA' || data.type === 'RNA') {
                    // Sadece DNA/RNA için ORF sonuçlarını göster
                    if (data.orfs && data.orfs.length > 0) {
                        // Kullanılan min ORF uzunluğunu başlığa yaz
                        minOrfLengthSpan.textContent = data.analysis_options_used?.min_orf_length || '?';
                        data.orfs.forEach(orf => {
                            const row = orfTableBody.insertRow();
                            row.classList.add('hover:bg-indigo-50'); // Satır üzerine gelince hafif vurgu
                            // Proteini kısalt ama tamamını title'da göster
                            const proteinSeq = orf.protein_sequence || '';
                            const proteinSeqShort = proteinSeq.length > 40
                                ? proteinSeq.substring(0, 18) + '...' + proteinSeq.substring(proteinSeq.length - 18)
                                : proteinSeq;

                            // Pozisyonları 1-tabanlı göster
                            const startNtDisplay = orf.position !== undefined ? orf.position + 1 : 'N/A'; // Updated to use 'position' from annotation object
                            const endNtDisplay = orf.end_position !== undefined ? orf.end_position + 1 : 'N/A'; // Updated to use 'end_position'

                            row.innerHTML = `
                                 <td class="px-3 py-1.5 whitespace-nowrap">${orf.frame || 'N/A'}</td>
                                 <td class="px-3 py-1.5 whitespace-nowrap">${orf.strand || 'N/A'}</td>
                                 <td class="px-3 py-1.5 whitespace-nowrap font-mono">${startNtDisplay}</td>
                                 <td class="px-3 py-1.5 whitespace-nowrap font-mono">${endNtDisplay}</td>
                                 <td class="px-3 py-1.5 whitespace-nowrap text-center">${orf.length_aa || 'N/A'}</td>
                                 <td class="px-3 py-1.5 font-mono text-[0.7rem] leading-tight max-w-xs truncate" title="${proteinSeq}">${proteinSeqShort || 'N/A'}</td>
                             `;
                        });
                        orfResultsDiv.classList.remove('hidden');
                        noOrfFoundP.classList.add('hidden');
                    } else {
                        orfResultsDiv.classList.remove('hidden'); // Bölümü göster
                        noOrfFoundP.classList.remove('hidden'); // Mesajı göster
                        minOrfLengthSpan.textContent = data.analysis_options_used?.min_orf_length || '?'; // Başlığı güncelle
                    }
                } else {
                    // Protein veya Unknown ise ORF bölümünü tamamen gizle
                    orfResultsDiv.classList.add('hidden');
                    noOrfFoundP.classList.add('hidden'); // Mesajı da gizle
                }
            }

            // --- 5. Protein İstatistiklerini Göster (Eğer varsa) ---
            if (proteinStatsResultsDiv) {
                let statsToDisplay = null;
                if (data.type === 'Protein') {
                    statsToDisplay = data.protein_stats; // Doğrudan protein analizi
                } else if ((data.type === 'DNA' || data.type === 'RNA') && data.orfs && data.orfs.length > 0) {
                    // En uzun ORF'un istatistikleri (backend'den gelmeli)
                    statsToDisplay = data.protein_stats; // Backend 'protein_stats'ı en uzun ORF için doldurmalı
                }

                if (statsToDisplay) {
                    displayProteinStats(statsToDisplay);
                    // Başlığı güncelle
                    const protStatsTitle = proteinStatsResultsDiv.querySelector('h3');
                    if (protStatsTitle) {
                        if (data.type === 'Protein') {
                            protStatsTitle.textContent = 'Protein İstatistikleri';
                        } else {
                            protStatsTitle.textContent = 'Protein İstatistikleri (En Uzun ORF İçin)';
                        }
                    }
                } else {
                    proteinStatsResultsDiv.classList.add('hidden'); // Veri yoksa bölümü gizle
                }
            }


            // --- 6. Diğer Detaylı Analizleri Doldur ---
            let hasOtherResultsSectionContent = false; // Bu detaylardan en az biri varsa otherResultsDiv'i göster
            if (otherResultsDiv) {

                // Motifler
                if (motifResultsDiv && motifListUl && noMotifFoundP) {
                    motifListUl.innerHTML = ''; // Temizle
                    if (data.type !== 'Protein' && data.motifs && data.motifs.length > 0) {
                        data.motifs.forEach(motif => {
                            const li = document.createElement('li');
                            // Position ve end_position artık 0-tabanlı geliyor
                            const startPosDisplay = motif.position !== undefined ? motif.position + 1 : 'N/A';
                            const endPosDisplay = motif.end_position !== undefined && motif.end_position !== motif.position ? `-${motif.end_position + 1}` : '';

                            li.innerHTML = `<span class="font-semibold">${motif.type}</span>: Poz ${startPosDisplay}${endPosDisplay}, Dizi: <span class="font-mono bg-gray-100 px-1 rounded">${motif.sequence || 'N/A'}</span>`;
                            li.title = motif.desc || motif.type; // Açıklamayı title'a ekle
                            li.classList.add('cursor-help'); // Yardım imleci
                            motifListUl.appendChild(li);
                        });
                        motifResultsDiv.classList.remove('hidden');
                        noMotifFoundP.classList.add('hidden');
                        hasOtherResultsSectionContent = true;
                    } else {
                        motifResultsDiv.classList.add('hidden');
                        if (data.type === 'DNA' || data.type === 'RNA') { // DNA/RNA ise ve motif yoksa mesajı göster
                            noMotifFoundP.classList.remove('hidden'); // Mesajı göster
                            motifResultsDiv.classList.remove('hidden'); // Bölümü aç (mesaj için)
                            hasOtherResultsSectionContent = true; // Mesaj olsa bile bölümü göster
                        }
                    }
                }

                // Restriksiyon Alanları
                if (restrictionResultsDiv && restrictionTableContainer && noRestrictionFoundP) {
                    restrictionTableContainer.innerHTML = ''; // Temizle
                    const sites = data.restriction_sites;
                    if (data.type === 'DNA' && sites && typeof sites === 'object' && !sites.error && Object.keys(sites).length > 0) {
                        let tableHTML = '<table class="min-w-full divide-y divide-gray-200"><thead class="bg-gray-100 sticky top-0"><tr><th class="px-2 py-1 text-left font-medium text-gray-600 uppercase tracking-wider">Enzim</th><th class="px-2 py-1 text-left font-medium text-gray-600 uppercase tracking-wider">Pozisyonlar (1-tabanlı)</th></tr></thead><tbody class="bg-white divide-y divide-gray-200">';
                        // Enzimleri alfabetik sırala
                        const sortedEnzymes = Object.keys(sites).sort();
                        sortedEnzymes.forEach(enzyme => {
                            // Pozisyonlar zaten 0-tabanlı geliyor, 1-tabanlı göstermek için +1 ekle
                            const positions = sites[enzyme].map(p => p + 1).join(', ');
                            tableHTML += `<tr class="hover:bg-indigo-50"><td class="px-2 py-1 whitespace-nowrap font-semibold">${enzyme}</td><td class="px-2 py-1 whitespace-nowrap font-mono">${positions}</td></tr>`;
                        });
                        tableHTML += '</tbody></table>';
                        restrictionTableContainer.innerHTML = tableHTML;
                        restrictionResultsDiv.classList.remove('hidden');
                        noRestrictionFoundP.classList.add('hidden');
                        hasOtherResultsSectionContent = true;
                    } else {
                        restrictionResultsDiv.classList.add('hidden');
                        if (data.type === 'DNA') { // Sadece DNA'da anlamlı
                            noRestrictionFoundP.classList.remove('hidden'); // Mesajı göster
                            if (sites && sites.error) {
                                noRestrictionFoundP.textContent = `Hesaplanamadı: ${sites.error}`;
                            } else {
                                noRestrictionFoundP.textContent = 'Tanımlı enzimler için kesim alanı bulunamadı.';
                            }
                            restrictionResultsDiv.classList.remove('hidden'); // Bölümü aç
                            hasOtherResultsSectionContent = true; // Mesajı göster
                        }
                    }
                }

                // CpG Adaları
                if (cpgResultsDiv && cpgListUl && noCpgFoundP && cpgAlgorithmUsedSpan) {
                    cpgListUl.innerHTML = ''; // Temizle
                    if (data.type === 'DNA' && data.cpg_islands && data.cpg_islands.length > 0) {
                        // Kullanılan algoritmayı başlığa yaz
                        cpgAlgorithmUsedSpan.textContent = data.analysis_options_used?.cpg_algorithm || 'varsayılan';
                        data.cpg_islands.forEach(island => {
                            // Position ve end_position artık 0-tabanlı geliyor
                            const startPosDisplay = island.position !== undefined ? island.position + 1 : 'N/A';
                            const endPosDisplay = island.end_position !== undefined ? island.end_position + 1 : 'N/A';

                            const li = document.createElement('li');
                            li.innerHTML = `Poz: <span class="font-mono">${startPosDisplay} - ${endPosDisplay}</span> (Uz: ${island.length} bp), %GC: ${island.gc_content}, O/E: ${island.cpg_oe?.toFixed(2)}`;
                            li.title = island.desc || 'CpG Adası'; // Açıklama
                            li.classList.add('cursor-help');
                            cpgListUl.appendChild(li);
                        });
                        cpgResultsDiv.classList.remove('hidden');
                        noCpgFoundP.classList.add('hidden');
                        hasOtherResultsSectionContent = true;
                    } else {
                        cpgResultsDiv.classList.add('hidden');
                        if (data.type === 'DNA') { // Sadece DNA'da anlamlı
                            cpgAlgorithmUsedSpan.textContent = data.analysis_options_used?.cpg_algorithm || 'varsayılan';
                            noCpgFoundP.classList.remove('hidden'); // Mesajı göster
                            cpgResultsDiv.classList.remove('hidden'); // Bölümü aç
                            hasOtherResultsSectionContent = true; // Mesajı göster
                        }
                    }
                }

                // Tüm Okuma Çerçeveleri (Çeviriler)
                if (framesResultsDiv && framesContentDiv) {
                    framesContentDiv.innerHTML = ''; // Temizle
                    if (data.type !== 'Protein' && data.frames && typeof data.frames === 'object' && !data.frames.error && Object.keys(data.frames).length > 0) {
                        let framesHTML = '';
                        const frameOrder = ['+1', '+2', '+3', '-1', '-2', '-3']; // Standart sıra
                        frameOrder.forEach(frame => {
                            if (data.frames[frame] !== undefined) {
                                const frameSeq = data.frames[frame] || '(Boş)'; // Boşsa belirt
                                const maxLength = 120; // Gösterilecek maksimum uzunluk
                                let truncatedSeq;
                                if (frameSeq.length > maxLength) {
                                    const firstPart = frameSeq.substring(0, 55);
                                    const lastPart = frameSeq.substring(frameSeq.length - 55);
                                    truncatedSeq = `${firstPart}...[${frameSeq.length - 110} aa]...${lastPart}`;
                                } else {
                                    truncatedSeq = frameSeq;
                                }
                                framesHTML += `<div class="mb-1"><span class="font-semibold text-gray-700">${frame}:</span> <span class="block md:inline break-all">${truncatedSeq}</span></div>`;
                            }
                        });
                        framesContentDiv.innerHTML = framesHTML;
                        framesResultsDiv.classList.remove('hidden');
                        hasOtherResultsSectionContent = true;
                    } else {
                        framesResultsDiv.classList.add('hidden');
                        if (data.frames && data.frames.error) {
                            framesContentDiv.innerHTML = `<p class="text-red-500 italic text-sm">Okuma çerçeveleri alınamadı: ${data.frames.error}</p>`;
                            framesResultsDiv.classList.remove('hidden'); // Bölümü aç (hata için)
                            hasOtherResultsSectionContent = true;
                        } else if (data.type === 'DNA' || data.type === 'RNA') {
                            // Hata yok ama veri de yok (beklenmedik durum?)
                            framesContentDiv.innerHTML = `<p class="text-gray-500 italic text-sm">Okuma çerçevesi verisi bulunamadı.</p>`;
                            framesResultsDiv.classList.remove('hidden');
                            hasOtherResultsSectionContent = true;
                        }
                    }
                }

                // Eğer detaylı analizlerden en least biri varsa veya mesaj gösteriliyorsa, ana bölümü göster
                if (hasOtherResultsSectionContent) {
                    otherResultsDiv.classList.remove('hidden');
                    // Accordion'ı açık başlat (isteğe bağlı)
                    const detailsElement = otherResultsDiv.querySelector('details');
                    if (detailsElement) detailsElement.open = true;
                } else {
                    otherResultsDiv.classList.add('hidden');
                }
            } // otherResultsDiv check sonu


            // --- Sonuç Container'ını Görünür Yap ---
            // Kısa bir gecikme ile opacity transition'ı tetikle
            requestAnimationFrame(() => {
                setTimeout(() => {
                    if (resultsContainer) resultsContainer.classList.add('visible-now');
                }, 10); // Çok kısa bir gecikme
            });


        } catch (error) {
            console.error("Sonuç Görüntüleme Hatası:", error);
            showError(`Sonuçlar işlenirken bir hata oluştu: ${error.message}. Konsolu kontrol edin.`);
            // Hata durumunda bile container'ı görünür yapabiliriz ki hata mesajı görünsün
            if (resultsContainer) resultsContainer.classList.add('visible-now');
        }
    } // displayResults sonu


    // --- Olay Dinleyicileri (Event Listeners) ---

    // Analiz Formu Gönderimi
    if (analysisForm) {
        analysisForm.addEventListener('submit', async (event) => {
            event.preventDefault();
            if (!sequenceInput || !analyzeButton || !loadingIndicator) {
                console.error("Form elemanları eksik!");
                showError("Arayüz hatası: Form elemanları bulunamadı.");
                return;
            }
            const sequenceText = sequenceInput.value.trim();

            if (!sequenceText) {
                showError("Lütfen analiz edilecek bir sekans girin veya yükleyin.");
                sequenceInput.focus();
                return;
            }

            showLoading();

            // Ayarları al
            const minOrfLength = orfLengthSlider ? parseInt(orfLengthSlider.value, 10) : 50;
            const cpgAlgorithm = cpgAlgorithmSelect ? cpgAlgorithmSelect.value : 'gardiner';

            const options = {
                min_orf_length: minOrfLength,
                cpg_algorithm: cpgAlgorithm
                // Gelecekte buraya başka analiz seçenekleri eklenebilir (örn. restriksiyon enzimleri, motif listesi vb.)
            };

            const payload = {
                sequence: sequenceText,
                options: options // Seçenekleri payload'a ekle
            };

            try {
                const response = await fetch('/analyze-sequence', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify(payload),
                });

                const data = await response.json();

                if (!response.ok) {
                    // Backend'den gelen hata mesajını kullan
                    throw new Error(data.error || `Sunucu hatası: ${response.status} ${response.statusText}`);
                }

                // Yanıtta backend hatası var mı diye kontrol et (örn. Python içinde yakalanan hata)
                if (data.error) {
                    throw new Error(`Analiz hatası: ${data.error}`);
                }

                // Her şey yolundaysa sonuçları göster
                displayResults(data);

            } catch (error) {
                console.error("Analiz İsteği veya İşleme Hatası:", error);
                // Kullanıcıya daha anlaşılır bir mesaj göster
                showError(`Analiz başarısız oldu. ${error.message}`);

            } finally {
                hideLoading(); // Her durumda yükleniyor göstergesini kaldır
            }
        });
    } else {
        console.error("Analysis form not found!");
    }

    // Dosya Yükleme
    if (fileInput && sequenceInput && fileNameSpan) {
        fileInput.addEventListener('change', (event) => {
            const file = event.target.files[0];
            if (file) {
                const reader = new FileReader();
                reader.onload = (e) => {
                    sequenceInput.value = e.target.result;
                    fileNameSpan.textContent = file.name;
                    fileNameSpan.title = file.name; // Tam adı tooltip olarak ekle
                    // Dosya yüklendiğinde otomatik analiz tetiklensin mi? İsteğe bağlı.
                    // analysisForm.dispatchEvent(new Event('submit'));
                };
                reader.onerror = (e) => {
                    showError(`Dosya okuma hatası: ${e.target.error}`);
                    fileNameSpan.textContent = '';
                    fileNameSpan.title = '';
                };
                reader.readAsText(file);
            } else {
                fileNameSpan.textContent = '';
                fileNameSpan.title = '';
            }
            // Input değerini temizle ki aynı dosya tekrar seçilebilsin
            event.target.value = null;
        });
    }

    // Örnek Sekans Butonları
    if (sampleButtons.length > 0 && sequenceInput && fileNameSpan) {
        sampleButtons.forEach(button => {
            button.addEventListener('click', () => {
                const seqId = button.dataset.seqId;
                if (sampleSequences[seqId]) {
                    sequenceInput.value = sampleSequences[seqId];
                    fileNameSpan.textContent = `Örnek: ${button.textContent}`; // Dosya adı yerine örnek adı
                    fileNameSpan.title = '';
                    // Örnek yüklendiğinde otomatik analiz?
                    // analysisForm.dispatchEvent(new Event('submit'));
                }
            });
        });
    }

    // Görselleştirme Alanına Tıklama/Hover (Event Delegation)
    // Mouseover ve mouseleave olayları renderVisualization fonksiyonunun içinde atanıyor.
    // Burada sadece detay gösterme/gizleme fonksiyonlarını tanımlayalım.
    function showAnnotationDetails(annotation) {
        if (!annotationDetailsDiv || !annotationDetailsContainer) return;

        // Detayları doldur
        let detailHtml = `<strong class="text-indigo-800">${annotation.type}</strong><br>`;

        // Pozisyonları 1-tabanlı göster
        const startPosDisplay = annotation.position !== undefined ? annotation.position + 1 : 'N/A';
        const endPosDisplay = annotation.end_position !== undefined && annotation.end_position !== annotation.position ? ` - ${annotation.end_position + 1}` : '';
        detailHtml += `Poz: ${startPosDisplay}${endPosDisplay}`;

        // Uzunluk
        if (annotation.length_aa) detailHtml += ` (Uz: ${annotation.length_aa} aa)`;
        else if (annotation.length) detailHtml += ` (Uz: ${annotation.length} bp)`;
        else if (annotation.end_position !== undefined) { // length yoksa hesapla
            // End position dahil olduğu için +1
            const len = Math.abs(annotation.end_position - annotation.position) + 1;
            if (len > 1) detailHtml += ` (Uz: ${len} ${currentSequenceType === 'Protein' ? 'aa' : 'bp'})`;
        }

        // Diğer bilgiler
        if (annotation.enzyme) detailHtml += `, Enzim: ${annotation.enzyme}`;
        if (annotation.strand) detailHtml += `, Strand: ${annotation.strand}`;
        if (annotation.frame) detailHtml += `, Çerçeve: ${annotation.frame}`;
        if (annotation.algorithm) detailHtml += `, Algoritma: ${annotation.algorithm}`;
        if (annotation.gc !== undefined) detailHtml += `, %GC: ${annotation.gc}`; // CpG Adası için GC
        if (annotation.oe !== undefined) detailHtml += `, O/E: ${annotation.oe?.toFixed(2)}`; // CpG Adası için O/E


        // Sekans (kısa ise)
        if (annotation.sequence) {
            const seqToDisplay = annotation.sequence.length < 50 ? annotation.sequence : annotation.sequence.substring(0, 45) + '...';
            detailHtml += `<br>Dizi: <span class="font-mono text-xs bg-white px-1 rounded">${seqToDisplay}</span>`;
        }


        // Eğitimsel Açıklama
        if (annotation.desc) {
            detailHtml += `<br><hr class="my-1 border-indigo-100"><span class="text-xs text-indigo-700">${annotation.desc}</span>`;
        }

        annotationDetailsDiv.innerHTML = detailHtml;
        annotationDetailsContainer.classList.remove('hidden');
    }

    function hideAnnotationDetails() {
        if (annotationDetailsContainer) {
            annotationDetailsContainer.classList.add('hidden');
        }
    }

    // Sayfa yüklendiğinde tema ve slider gibi başlangıç ayarları
    // Bu kısım index.html içinden script.js'in sonuna taşındı veya kopyalandı.
    // İki kez çalışmasını önlemek için index.html'deki kopyayı kaldırın.
    // veya burada tanımlanan fonksiyonları dışa aktarıp index.html'den çağırın.
    // Şimdilik index.html'deki script bloğunun kaldırıldığını varsayalım.

    const themeToggle = document.getElementById('theme-toggle');
    const body = document.getElementById('main-body'); // index.html'deki body'ye id eklenmeli

    // Kaydedilmiş tema tercihini kontrol et ve uygula
    const savedTheme = localStorage.getItem('theme');
    const prefersDark = window.matchMedia('(prefers-color-scheme: light)').matches;
    let isDarkMode;

    if (savedTheme) {
        isDarkMode = savedTheme === 'light';
    } else {
        isDarkMode = prefersDark; // Sistem tercihini kullan
        //localStorage.setItem('theme', isDarkMode ? 'dark' : 'light'); // İlk yüklemede kaydetmeyelim, kullanıcı seçerse kaydederiz
    }

    if (isDarkMode) {
        body.classList.add('dark-mode');
    }
    updateThemeIcon(isDarkMode);


    // Tema değiştirme butonu işlevi
    if (themeToggle && body) {
        themeToggle.addEventListener('click', () => {
            const newIsDarkMode = body.classList.toggle('dark-mode');
            localStorage.setItem('theme', newIsDarkMode ? 'dark' : 'light');
            updateThemeIcon(newIsDarkMode);
        });
    }


    function updateThemeIcon(isDarkMode) {
        if (!themeToggle) return;
        // İçerideki SVG'yi tamamen değiştirmek yerine, class ekleyip CSS ile görünürlüğü kontrol edebiliriz
        // veya doğrudan SVG path'ini değiştirebiliriz. Şimdilik SVG'yi değiştirelim.
        if (isDarkMode) {
            // Sun icon (açık tema ikonu)
            themeToggle.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" class="h-6 w-6" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2">
                <path stroke-linecap="round" stroke-linejoin="round" d="M12 3v1m0 16v1m9-9h-1M4 12H3m15.364 6.364l-.707-.707M6.343 6.343l-.707-.707m12.728 0l-.707.707M6.343 17.657l-.707.707M16 12a4 4 0 11-8 0 4 4 0 018 0z" />
            </svg>`;
        } else {
            // Moon icon (koyu tema ikonu)
            themeToggle.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" class="h-6 w-6" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2">
                <path stroke-linecap="round" stroke-linejoin="round" d="M20.354 15.354A9 9 0 018.646 3.646 9.003 9.003 0 0012 21a9.003 9.003 0 008.354-5.646z" />
            </svg>`;
        }
    }

    // ORF uzunluk slider değeri gösterimi
    const orfLengthSliderEl = document.getElementById('orf-length-slider');
    const orfLengthValueEl = document.getElementById('orf-length-value');

    if (orfLengthSliderEl && orfLengthValueEl) {
        // Initial value update
        orfLengthValueEl.textContent = `${orfLengthSliderEl.value} aa`;
        // Update on input
        orfLengthSliderEl.addEventListener('input', () => {
            orfLengthValueEl.textContent = `${orfLengthSliderEl.value} aa`;
        });
    }

    // Global Chart objesini tutmak için değişken
    window.distributionChart = null;

    // Dağılım grafiği oluşturma fonksiyonu
    window.createDistributionChart = function (distribution, seqType) {
        const ctx = document.getElementById('distribution-chart');
        // Ensure the parent container is ready to hold the canvas
        const chartContainer = ctx ? ctx.parentElement : null;
        if (!ctx || !chartContainer || !distribution || Object.keys(distribution).length === 0) {
            // Eğer canvas veya veri yoksa veya önceki chart varsa temizle
            if (window.distributionChart) {
                window.distributionChart.destroy();
                window.distributionChart = null;
            }
            // Grafik yerine mesaj gösterilebilir
            if (chartContainer) {
                chartContainer.innerHTML = '<p class="text-center text-xs text-gray-500 italic p-4">Dağılım verisi yok veya yetersiz.</p>';
                // Restore the canvas element if it was removed
                if (!chartContainer.querySelector('canvas#distribution-chart')) {
                    const newCanvas = document.createElement('canvas');
                    newCanvas.id = 'distribution-chart';
                    newCanvas.height = 256; // Match container height or set standard
                    chartContainer.appendChild(newCanvas);
                }
            }
            return;
        }

        // Önceki grafik varsa yok et
        if (window.distributionChart) {
            window.distributionChart.destroy();
        }

        const labels = Object.keys(distribution).sort(); // Alfabetik sırala
        const data = labels.map(label => distribution[label]); // Sıralanmış etiketlere göre veriyi al

        // Renkleri belirle (seqType'a göre)
        let backgroundColors = [];
        let borderColors = [];
        const dnaRnaColorMap = {
            'A': 'rgba(52, 211, 153, 0.7)', // Emerald-400
            'T': 'rgba(248, 113, 113, 0.7)', // Red-400
            'U': 'rgba(248, 113, 113, 0.7)', // Red-400 (RNA)
            'G': 'rgba(251, 191, 36, 0.7)',  // Amber-400
            'C': 'rgba(96, 165, 250, 0.7)',  // Blue-400
        };
        const proteinColorMap = {
            // Basit bir harita veya gruplama kullanılabilir
            'A': '#3498db', 'V': '#2980b9', 'L': '#1f618d', 'I': '#1a5276', // Hydrophobic shades
            'M': '#154360', 'F': '#0e2f44', 'W': '#0a2229', 'P': '#34495e', 'G': '#7f8c8d', // Special/Other Hydrophobic
            'S': '#2ecc71', 'T': '#27ae60', 'N': '#229954', 'Q': '#1e8449', // Polar shades
            'K': '#e74c3c', 'R': '#c0392b', 'H': '#a93226', // Positive shades
            'D': '#f39c12', 'E': '#d35400', // Negative shades
            'C': '#95a5a6', // Special
            'X': '#cccccc', '*': '#999999', '-': '#666666', // Unknown/Stop/Gap
            'default': 'rgba(156, 163, 175, 0.7)' // Default gray
        };

        if (seqType === 'DNA' || seqType === 'RNA') {
            backgroundColors = labels.map(label => dnaRnaColorMap[label] || dnaRnaColorMap['default'] || 'rgba(156, 163, 175, 0.7)');
        } else if (seqType === 'Protein') {
            backgroundColors = labels.map(label => proteinColorMap[label] || proteinColorMap['default'] || 'rgba(156, 163, 175, 0.7)');
        } else {
            backgroundColors = labels.map(label => 'rgba(156, 163, 175, 0.7)'); // Generic gray for unknown
        }

        borderColors = backgroundColors.map(color => {
            if (color.startsWith('rgba')) {
                return color.replace('0.7', '1'); // Make border opaque
            }
            return color; // Assume solid color if not rgba
        });


        window.distributionChart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: labels,
                datasets: [{
                    label: 'Dağılım (%)',
                    data: data,
                    backgroundColor: backgroundColors,
                    borderColor: borderColors,
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false, // Container'a uyması için önemli
                indexAxis: 'x', // Yatay çubuklar için 'y'
                scales: {
                    y: {
                        beginAtZero: true,
                        title: { display: true, text: 'Yüzde (%)' },
                        grid: { color: getComputedStyle(document.body).getPropertyValue('--grid-color') || 'rgba(200, 200, 200, 0.1)' } // CSS değişkeninden al
                    },
                    x: {
                        title: { display: true, text: (seqType === 'Protein' ? 'Amino Asit' : 'Baz') },
                        grid: { display: false } // X ekseni gridini kapat
                    }
                },
                plugins: {
                    legend: {
                        display: false // Tek dataset olduğu için gereksiz
                    },
                    tooltip: {
                        backgroundColor: getComputedStyle(document.body).getPropertyValue('--tooltip-bg-color') || 'rgba(0, 0, 0, 0.7)',
                        titleFont: { size: 14 },
                        bodyFont: { size: 12 },
                        callbacks: {
                            label: function (context) {
                                let label = context.dataset.label || '';
                                if (label) { label += ': '; }
                                if (context.parsed.y !== null) {
                                    // Yüzde formatı
                                    label += `${context.parsed.y.toFixed(2)}%`;
                                }
                                return label;
                            }
                        }
                    }
                },
                animation: {
                    duration: 500 // Hafif animasyon
                }
            }
        });
    }; // createDistributionChart sonu


}); // DOMContentLoaded sonu
