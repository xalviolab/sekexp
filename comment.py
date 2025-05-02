# -*- coding: utf-8 -*-
"""
Bu modül, Sekans Explorer uygulamasında gösterilen çeşitli biyoinformatik analiz
sonuçları için eğitimsel açıklamalar ve bağlam sağlar.

Amacı, kullanıcılara analiz edilen özelliklerin (örn. GC içeriği, kodonlar, motifler)
ne anlama geldiğini, nasıl hesaplandığını ve potansiyel biyolojik önemini
anlaşılır bir dilde sunmaktır.

Bu açıklamalar genel bilgi amaçlıdır ve spesifik araştırma veya klinik
kararlar için temel alınmamalıdır.
"""

# Açıklamalar genellikle bir sözlük yapısında saklanır.
# Anahtar: Analiz edilen özelliğin kodu (örn. 'gc_content', 'start_codon')
# Değer: Açıklama metni veya daha yapılandırılmış bir bilgi (örn. başlık, detay, referans)
EDUCATIONAL_COMMENTS = {
    # --- Genel İstatistikler ---
    'sequence_length': {
        'title': "Sekans Uzunluğu",
        'explanation': "Analiz edilen DNA, RNA veya protein dizisinin toplam nükleotid veya amino asit sayısıdır.",
        'calculation': "Dizideki karakterlerin sayılmasıyla belirlenir.",
        'context': "Uzunluk, genin veya proteinin büyüklüğü hakkında genel bir fikir verir."
    },
    'gc_content': {
        'title': "GC İçeriği (%)",
        'explanation': "DNA veya RNA dizisindeki Guanin (G) ve Sitozin (C) bazlarının toplam oranını yüzde olarak ifade eder.",
        'calculation': "(G sayısı + C sayısı) / Toplam baz sayısı * 100 formülü ile hesaplanır.",
        'context': "Yüksek GC içeriği, DNA'nın erime sıcaklığını (Tm) artırabilir ve gen regülasyonu bölgeleriyle (CpG adaları gibi) ilişkili olabilir."
    },
    'at_content': {
        'title': "AT İçeriği (%)",
        'explanation': "DNA veya RNA dizisindeki Adenin (A) ve Timin (T) veya Urasil (U) bazlarının toplam oranını yüzde olarak ifade eder.",
        'calculation': "(A sayısı + T/U sayısı) / Toplam baz sayısı * 100 veya 100 - GC İçeriği (%) formülü ile hesaplanır.",
        'context': "Genellikle GC içeriği ile ters orantılıdır. Bazı DNA bölgeleri (örn. TATA kutusu) AT bakımından zengindir."
    },
    'molecular_weight_dna_rna': {
        'title': "Molekül Ağırlığı (DNA/RNA)",
        'explanation': "DNA veya RNA dizisinin tahmini molekül ağırlığıdır (genellikle g/mol cinsinden).",
        'calculation': "Dizideki her bir nükleotidin ortalama molekül ağırlığına göre hesaplanır (çift sarmal veya tek sarmal varsayımına göre değişebilir).",
        'context': "Molekül ağırlığı, laboratuvar tekniklerinde (örn. jel elektroforezi) ve kantitatif analizlerde kullanılır."
    },
    'molecular_weight_protein': {
        'title': "Molekül Ağırlığı (Protein)",
        'explanation': "Protein dizisinin tahmini molekül ağırlığıdır (genellikle Dalton (Da) veya kiloDalton (kDa) cinsinden).",
        'calculation': "Dizideki her bir amino asidin ortalama molekül ağırlığının toplanmasıyla hesaplanır.",
        'context': "Protein büyüklüğünü gösterir ve saflaştırma, tanımlama (örn. Western Blot) gibi tekniklerde önemlidir."
    },
    'tm_estimate': {
        'title': "Ergime Sıcaklığı (Tm) Tahmini",
        'explanation': "DNA çift sarmalının yarısının tek sarmallara ayrıldığı tahmini sıcaklıktır (°C).",
        'calculation': "Basit formüllerle (Wallace kuralı, GC içeriği tabanlı) veya daha karmaşık modellerle tahmin edilir. Bu analizde kullanılan yöntem basit bir yaklaşımdır.",
        'context': "PCR primer tasarımı, hibridizasyon deneyleri gibi moleküler biyoloji uygulamalarında önemlidir. Tahminler, tuz konsantrasyonu gibi faktörlerden etkilenebilir."
    },

    # --- DNA/RNA Özellikleri ---
    'start_codon': {
        'title': "Başlangıç Kodonu (örn. ATG)",
        'explanation': "Protein sentezinin (translasyon) başladığı üç nükleotidlik kodondur. Genellikle Metiyonin amino asidini kodlar.",
        'calculation': "Dizide 'ATG' (DNA'da) veya 'AUG' (RNA'da) dizilerinin aranmasıyla bulunur.",
        'context': "Açık Okuma Çerçevesi'nin (ORF) başlangıcını işaret eder. Ancak her ATG/AUG kodonu translasyonu başlatmayabilir, bağlam önemlidir."
    },
    'stop_codon': {
        'title': "Bitiş Kodonu (örn. TAA, TAG, TGA)",
        'explanation': "Protein sentezinin sonlandığı üç nükleotidlik kodonlardır. Bu kodonlar bir amino asit kodlamaz.",
        'calculation': "Dizide 'TAA', 'TAG', 'TGA' (DNA'da) veya 'UAA', 'UAG', 'UGA' (RNA'da) dizilerinin aranmasıyla bulunur.",
        'context': "Açık Okuma Çerçevesi'nin (ORF) sonunu belirler."
    },
    'TATA-box': {
        'title': "TATA Kutusu",
        'explanation': "Ökaryotik genlerin promotor bölgesinde bulunan, transkripsiyonun başlamasında kritik rol oynayan konsensüs bir DNA dizisidir (genellikle TATAAA).",
        'calculation': "Dizide TATA[AT]A[AT] gibi motiflerin aranmasıyla bulunur.",
        'context': "RNA polimeraz II ve diğer transkripsiyon faktörlerinin bağlanması için bir platform görevi görür."
    },
    'CpG': {
        'title': "CpG Dinükleotidi",
        'explanation': "Bir Sitozin (C) bazının hemen ardından bir Guanin (G) bazının geldiği ('C' fosfat 'G') dizidir.",
        'calculation': "Dizide 'CG' ikililerinin sayılmasıyla bulunur.",
        'context': "Memeli genomlarında CpG dinükleotidleri genellikle metilasyona uğrar. CpG adaları (CpG'lerin yoğun olduğu bölgeler) genellikle gen promotorlarında bulunur ve metilasyon durumu gen ifadesini etkiler."
    },
    'cpg_island': {
        'title': "CpG Adası",
        'explanation': "Genomda CpG dinükleotidlerinin beklenenden daha yüksek sıklıkta bulunduğu bölgelerdir. Genellikle genlerin başlangıç bölgeleriyle (promotorlar) ilişkilidir.",
        'calculation': "Belirli bir uzunluktaki (örn. >200 bp) pencerede GC içeriğinin (>%50) ve gözlenen CpG oranının beklenene göre oranının (>0.6) belirli eşik değerlerini aşmasıyla tanımlanır.",
        'context': "CpG adaları genellikle metillenmemiş durumdadır ve aktif genlerle ilişkilidir. Metilasyonları gen susturulmasına yol açabilir."
    },
    'ORF': {
        'title': "Açık Okuma Çerçevesi (ORF)",
        'explanation': "Bir başlangıç kodonu (örn. ATG) ile başlayıp bir bitiş kodonu (örn. TAA, TAG, TGA) ile biten, potansiyel olarak bir protein kodlayabilen DNA veya RNA dizisi bölümüdür.",
        'calculation': "Altı olası okuma çerçevesinde (üçü ileri, üçü geri yönde) başlangıç ve bitiş kodonları arasındaki dizilerin belirlenmesiyle bulunur.",
        'context': "Gen bulma ve protein kodlama potansiyelini analiz etmede temel bir adımdır. Uzun ORF'ler genellikle fonksiyonel genleri işaret eder."
    },
    'restriction_site': {
        'title': "Restriksiyon Bölgesi",
        'explanation': "Restriksiyon enzimlerinin tanıdığı ve spesifik olarak kestiği kısa DNA dizileridir.",
        'calculation': "Bilinen restriksiyon enzimlerinin tanıma dizilerinin (örn. EcoRI için GAATTC) sekans içinde aranmasıyla bulunur.",
        'context': "Moleküler klonlama, gen haritalama ve RFLP analizi gibi birçok moleküler biyoloji tekniğinde kullanılır."
    },

    # --- Protein Özellikleri ---
    'isoelectric_point': {
        'title': "İzoelektrik Nokta (pI)",
        'explanation': "Bir proteinin net elektrik yükünün sıfır olduğu pH değeridir.",
        'calculation': "Proteindeki asidik ve bazik amino asit kalıntılarının pKa değerleri kullanılarak iteratif yöntemlerle hesaplanır.",
        'context': "Protein saflaştırma (izoelektrik odaklama gibi), çözünürlük ve protein-protein etkileşimlerini anlamada önemlidir."
    },
    'amino_acid_composition': {
        'title': "Amino Asit Kompozisyonu",
        'explanation': "Proteindeki her bir amino asit türünün sayısını ve yüzdesini gösterir.",
        'calculation': "Dizideki her bir amino asit karakterinin sayılmasıyla belirlenir.",
        'context': "Protein özellikleri (örn. hidrofobisite, yük), evrimsel analizler ve protein mühendisliği için bilgi sağlar."
    },
    'instability_index': {
        'title': "İnstabilite İndeksi",
        'explanation': "Bir proteinin in vitro (test tüpünde) kararlılığı hakkında bir tahmin sunar. 40'ın üzerindeki değerler proteinin kararsız (instabil) olabileceğini düşündürür.",
        'calculation': "Belirli dipeptitlerin protein dizisindeki ağırlıklı frekanslarına dayanarak ampirik bir formülle hesaplanır.",
        'context': "Rekombinant protein üretimi ve proteinin yarı ömrü hakkında fikir verebilir, ancak in vivo (canlı hücredeki) kararlılığı tam olarak yansıtmayabilir."
    },
    'gravy_index': {
        'title': "GRAVY İndeksi (Grand Average of Hydropathicity)",
        'explanation': "Protein dizisindeki tüm amino asitlerin hidropati (su sevme/sevmeme) değerlerinin ortalamasıdır.",
        'calculation': "Her amino asidin hidropati değeri toplanıp toplam amino asit sayısına bölünerek hesaplanır.",
        'context': "Pozitif değerler genellikle hidrofobik (suyu sevmeyen, örn. membran proteinleri), negatif değerler ise hidrofilik (suyu seven, örn. çözünür proteinler) proteinleri işaret eder."
    },
    'secondary_structure_prediction': {
        'title': "İkincil Yapı Tahmini",
        'explanation': "Protein dizisinin yerel olarak hangi ikincil yapı elemanlarını (alfa-heliks, beta-tabaka, koil/döngü) oluşturma eğiliminde olduğunu tahmin eder.",
        'calculation': "Genellikle amino asitlerin istatistiksel eğilimlerine veya makine öğrenmesi modellerine dayalı algoritmalarla (örn. Chou-Fasman, GOR, PSIPRED) yapılır. Bu analizdeki tahmin basit bir yaklaşımdır.",
        'context': "Proteinlerin üç boyutlu yapısını ve fonksiyonunu anlamada ilk adımdır. Tahminlerin doğruluğu değişkendir."
    },

    # --- Motifler ve Diğer Özellikler (sequence_analyzer.py'den alınabilir) ---
    # ANNOTATION_TEMPLATES içindeki diğer anahtarlar buraya eklenebilir.
    # Örnek:
    'promoter_-10': {
        'title': "Promotör -10 Kutusu (Pribnow)",
        'explanation': "Bakteriyel promotörlerde transkripsiyon başlangıç noktasının ~10 baz öncesinde bulunan, RNA polimerazın bağlanmasında rol oynayan TATAAT benzeri dizi.",
        'calculation': "Dizide TATAAT ve benzeri motiflerin aranmasıyla bulunur.",
        'context': "Bakteriyel gen ifadesinin düzenlenmesinde önemlidir."
    },
    # ... Diğer motifler ...

}

def get_comment(feature_key):
    """Belirli bir özellik anahtarı için eğitimsel yorumu döndürür."""
    return EDUCATIONAL_COMMENTS.get(feature_key, {
        'title': feature_key.replace('_', ' ').title(),
        'explanation': "Bu özellik için detaylı açıklama henüz mevcut değil.",
        'calculation': "-",
        'context': "-"
    })

# İleride eklenebilecekler:
# - Motiflerin hangi protein ailelerinde görüldüğüne dair bilgiler (örn. Pfam, InterPro bağlantıları)
# - İlgili literatür referansları
# - Daha interaktif açıklamalar için fonksiyonlar