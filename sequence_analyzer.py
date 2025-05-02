import re
from Bio.Seq import Seq
from collections import Counter
import numpy as np  # GC veya diğer hesaplamalar için potansiyel olarak kullanılabilir
from Bio.SeqUtils import gc_fraction as GC_biopython # Rename to avoid conflict with GC-box
from Bio.SeqUtils import ProtParam  # Protein analizleri için
from Bio.Restriction import Analysis, RestrictionBatch, CommOnly # Restriksiyon analizi için
import traceback

import app # For detailed error logging


# --- Eğitimsel Açıklama Şablonları ---
# DİKKAT: Bu açıklamalar tamamen eğitimseldir ve ASLA tıbbi/klinik anlam taşımaz.
ANNOTATION_TEMPLATES = {
    # Kodonlar
    "start_codon": "Başlangıç Kodonu ({position}): Genellikle Metiyonin (M) amino asidini kodlar ve protein sentezinin başladığı muhtemel noktayı işaret eder.",
    "stop_codon": "Bitiş Kodonu ({position}): Protein sentezinin sonlandığı sinyaldir. Bu kodonlar bir amino asit kodlamaz.",

    # Temel promotör elementleri
    "TATA-box": "TATA Kutusu ({position}): Dizi: '{sequence_context}' - Ökaryotik genlerin promotor bölgesinde bulunan, transkripsiyonun başlamasında rol oynayan yaygın bir DNA dizisidir.",
    "CAAT-box": "CAAT Kutusu ({position}): Dizi: '{sequence_context}' - Ökaryotik genlerin promotor bölgesinde bulunan, transkripsiyonun düzenlenmesinde rol oynayan yaygın bir DNA dizisidir.",
    "GC-box": "GC Kutusu ({position}): Dizi: '{sequence_context}' - Ökaryotik genlerin promotor bölgesinde bulunan, transkripsiyonun düzenlenmesinde rol oynayan GC bakımından zengin bir DNA dizisidir.",
    "BRE": "TFIIB Tanıma Elementi (BRE) ({position}): Dizi: '{sequence_context}' - TATA kutusunun yukarı akışında bulunan, TFIIB transkripsiyon faktörünün bağlandığı bölgedir.",
    "Inr": "Başlatıcı Element (Inr) ({position}): Dizi: '{sequence_context}' - Transkripsiyon başlangıç noktasında bulunan, RNA polimeraz II'nin bağlanmasında rol oynayan dizisidir.",
    "DPE": "Aşağı Akış Promotör Elementi (DPE) ({position}): Dizi: '{sequence_context}' - Transkripsiyon başlangıç noktasının yaklaşık 30 baz çifti aşağısında bulunan, TFIID bağlanmasında rol oynayan elementtir.",

    # Transkripsiyon düzenleyici elementler
    "enhancer_motif": "Güçlendirici Motif ({position}): Dizi: '{sequence_context}' - Gen ifadesinin seviyesini artırabilen düzenleyici bir DNA bölgesidir.",
    "silencer_motif": "Susturucu Motif ({position}): Dizi: '{sequence_context}' - Gen ifadesinin seviyesini azaltabilen veya baskılayabilen düzenleyici bir DNA bölgesidir.",
    "AP-1_site": "AP-1 Bağlanma Bölgesi ({position}): Dizi: '{sequence_context}' - AP-1 transkripsiyon faktörünün bağlandığı, gen ifadesinin düzenlenmesinde rol oynayan TGA[CG]TCA motifi.",
    "NF-kB_site": "NF-kB Bağlanma Bölgesi ({position}): Dizi: '{sequence_context}' - NF-kB transkripsiyon faktörünün bağlandığı, bağışıklık ve inflamasyon yanıtlarında önemli olan GGG[AG][CT]T[CT][CT]CC motifi.",
    "E-box": "E-Kutusu Motifi ({position}): Dizi: '{sequence_context}' - bHLH (basic Helix-Loop-Helix) transkripsiyon faktörlerinin bağlandığı CANNTG motifi (N=herhangi bir baz).",
    "HRE": "Hipoksi Yanıt Elementi (HRE) ({position}): Dizi: '{sequence_context}' - Düşük oksijen koşullarında HIF-1 transkripsiyon faktörünün bağlandığı bölgedir.",
    "SRE": "Serum Yanıt Elementi (SRE) ({position}): Dizi: '{sequence_context}' - Serum uyarısına yanıt olarak gen ifadesini düzenleyen, SRF transkripsiyon faktörünün bağlandığı bölgedir.",
    "CArG-box": "CArG Kutusu ({position}): Dizi: '{sequence_context}' - MADS-box transkripsiyon faktörlerinin bağlandığı, CC(A/T)6GG konsensüs dizisini içeren bölgedir.",
    "GATA_motif": "GATA Motifi ({position}): Dizi: '{sequence_context}' - GATA transkripsiyon faktörlerinin bağlandığı, özellikle kan hücrelerinin gelişiminde rol oynayan [AT]GATA[AG] motifi.",
    "CREB_site": "CREB Bağlanma Bölgesi (CRE) ({position}): Dizi: '{sequence_context}' - CREB transkripsiyon faktörünün bağlandığı, hücre sinyal iletiminde rol oynayan TGACGTCA motifi.",
    "heat_shock_element": "Isı Şoku Elementi (HSE) ({position}): Dizi: '{sequence_context}' - Isı şoku veya diğer stres koşullarına yanıt olarak aktive olan genlerin promotörlerinde bulunan GAANNTTC benzeri motif.",
    "STAT_binding_site": "STAT Bağlanma Bölgesi ({position}): Dizi: '{sequence_context}' - STAT transkripsiyon faktörlerinin bağlandığı, sitokin sinyallerinin aktarımında rol oynayan bölgedir.",
    "p53_binding_site": "p53 Bağlanma Bölgesi ({position}): Dizi: '{sequence_context}' - Tümör baskılayıcı p53 proteininin bağlandığı, hücre döngüsü kontrolü ve apoptoz ile ilişkili bölgedir.",
    "NFAT_binding_site": "NFAT Bağlanma Bölgesi ({position}): Dizi: '{sequence_context}' - NFAT transkripsiyon faktörünün bağlandığı, bağışıklık yanıtı ve T hücresi aktivasyonunda rol oynayan bölgedir.",
    "ISRE": "İnterferon-Uyarılı Yanıt Elementi (ISRE) ({position}): Dizi: '{sequence_context}' - İnterferon sinyallemesinde rol oynayan, IRF transkripsiyon faktörlerinin bağlandığı bölgedir.",

    # RNA işleme sinyalleri
    "poly_a_signal": "Poli-A Sinyali ({position}): Dizi: '{sequence_context}' - mRNA'nın 3' ucuna poli-A kuyruğunun eklenmesi için gerekli olan sinyal dizisidir.",
    "splice_donor": "Splice Donör (5' Splice Alanı) ({position}): Dizi: '{sequence_context}' - İntronların başlangıcında bulunan ve splaysing (kesip birleştirme) işlemi için önemli olan konsensüs dizisi.",
    "splice_acceptor": "Splice Akseptör (3' Splice Alanı) ({position}): Dizi: '{sequence_context}' - İntronların sonunda bulunan ve splaysing işlemi için önemli olan konsensüs dizisi.",
    "branch_point": "Dallanma Noktası ({position}): Dizi: '{sequence_context}' - İntronlarda bulunan, splicing işlemi sırasında lariat (kement) oluşumunda rol oynayan adenin içeren bölgedir.",

    # Bakteriyel elementler
    "Shine-Dalgarno": "Shine-Dalgarno Sırası ({position}): Dizi: '{sequence_context}' - Bakterilerde ribozomun mRNA'ya bağlanmasını sağlayan, translasyon başlangıcı öncesinde bulunan purin bakımından zengin dizidir.",
    "RBS": "Ribozom Bağlanma Bölgesi (RBS) ({position}): Dizi: '{sequence_context}' - Bakterilerde ribozomun mRNA'ya bağlanmasını sağlayan, translasyon başlangıcı öncesinde bulunan bölgedir.",
    "promoter_-10": "Promotör -10 Kutusu (Pribnow) ({position}): Dizi: '{sequence_context}' - Bakteriyel promotörlerde transkripsiyon başlangıç noktasının ~10 baz öncesinde bulunan, RNA polimerazın bağlanmasında rol oynayan TATAAT benzeri dizi.",
    "promoter_-35": "Promotör -35 Kutusu ({position}): Dizi: '{sequence_context}' - Bakteriyel promotörlerde transkripsiyon başlangıç noktasının ~35 baz öncesinde bulunan, RNA polimerazın tanınmasında rol oynayan TTGACA benzeri dizi.",
    "rho_independent_terminator": "Rho-Bağımsız Terminatör ({position}): Dizi: '{sequence_context}' - Bakteriyel transkripsiyon sonlandırma sinyali, GC bakımından zengin bir saç tokası yapısı ve ardından gelen T bakımından zengin bir bölgeden oluşur.",

    # DNA replikasyon elementleri
    "ori_consensus": "Replikasyon Orijini Konsensüsü ({position}): Dizi: '{sequence_context}' - DNA replikasyonunun başladığı bölgede bulunabilen tekrarlayan motif örneği.",
    "DnaA_box": "DnaA Kutusu ({position}): Dizi: '{sequence_context}' - Bakterilerde replikasyon orijinine DnaA proteininin bağlanmasını sağlayan TTATCCACA benzeri dizi.",
    "autonomously_replicating_sequence": "Otonom Replikasyon Dizisi (ARS) ({position}): Dizi: '{sequence_context}' - Mayada DNA replikasyonunun başlangıç noktası olarak işlev gören konsensüs dizisidir.",

    # Tekrar dizileri
    "simple_repeat_AT": "Basit AT Tekrarı ({position}): Dizi: '{sequence_context}' - AT dinükleotidinin ardışık tekrarlarından oluşan bölgedir.",
    "simple_repeat_GC": "Basit GC Tekrarı ({position}): Dizi: '{sequence_context}' - GC dinükleotidinin ardışık tekrarlarından oluşan bölgedir.",
    "simple_repeat_CA": "Basit CA Tekrarı ({position}): Dizi: '{sequence_context}' - CA dinükleotidinin ardışık tekrarlarından oluşan bölgedir.",
    "simple_repeat_GA": "Basit GA Tekrarı ({position}): Dizi: '{sequence_context}' - GA dinükleotidinin ardışık tekrarlarından oluşan bölgedir.",
    "trinucleotide_repeat": "Trinükleotid Tekrarı ({position}): Dizi: '{sequence_context}' - Üç nükleotidin ardışık tekrarlarından oluşan bölgedir, bazı nörodejeneratif hastalıklarla ilişkili olabilir.",

    # Telomer ve Sentromer Elementleri
    "telomeric_repeat": "Telomerik Tekrar ({position}): Dizi: '{sequence_context}' - Kromozomların uçlarında bulunan, TTAGGG dizisinin tekrarlarından oluşan bölgedir.",
    "CENP_B_box": "CENP-B Kutusu ({position}): Dizi: '{sequence_context}' - Sentromer protein B'nin bağlandığı, sentromer bölgesinde bulunan 17 bp'lik konsensüs dizidir.",

    # Metilasyon İlişkili Elementler
    "CpG": "CpG Dinükleotidi ({position}): Dizi: '{sequence_context}' - DNA metilasyonunun gerçekleştiği potansiyel bölge, memeli genomlarında nispeten az bulunur.",
    "CTCF_binding_site": "CTCF Bağlanma Bölgesi ({position}): Dizi: '{sequence_context}' - İzolatör proteini CTCF'nin bağlandığı, kromatin organizasyonunda ve gen ifadesinin düzenlenmesinde rol oynayan bölgedir.",

    # Mobil Genetik Elementler
    "LTR_consensus": "Uzun Terminal Tekrar (LTR) Konsensüsü ({position}): Dizi: '{sequence_context}' - Retrovirüsler ve retrotranspozonlarda bulunan, gen ifadesini düzenleyen tekrarlı dizilerdir.",
    "transposon_IR": "Transpozon Ters Tekrarı (IR) ({position}): Dizi: '{sequence_context}' - Transpozonların uçlarında bulunan, transpozaz enziminin tanıdığı ters tekrarlı dizilerdir.",

    # Diğer özellikler
    "cpg_island": "CpG Adası ({position}-{end_position}): Uz: {length} bp, %GC: {gc_content}, O/E: {cpg_oe:.2f} (Algoritma: {algorithm}) - Gen promotörlerinde sıklıkla bulunan, CG dinükleotidlerince zengin ve genellikle metilasyona uğramayan DNA bölgesidir.",
    "restriction_site": "Restriksiyon Bölgesi ({position}): Enzim: {enzyme} - DNA'yı spesifik olarak tanıyan ve kesen '{enzyme}' enziminin tanıma bölgesidir.",
    "ORF": "Açık Okuma Çerçevesi (ORF) ({position}-{end_position}): Uz: {length_aa} aa ({length_nt} bp), Çerçeve: {frame}, Strand: {strand} - Potansiyel protein kodlayan bölge, başlangıç kodonundan bitiş kodonuna kadar olan kısımdır."
}

# --- Yardımcı Fonksiyon: Eğitimsel Açıklama Üretimi ---
def generate_educational_annotation(feature_type, position, sequence_context=None, length=None, gc_content=None, cpg_oe=None, enzyme=None, length_aa=None, frame=None, end_position=None, length_nt=None, strand=None, algorithm=None, **kwargs):
    """Verilen özellik türü ve parametreler için eğitimsel bir açıklama üretir."""
    template = ANNOTATION_TEMPLATES.get(feature_type, "Tanımlanmamış Özellik ({position})")

    # Pozisyonları 1-tabanlı yapalım (biyologlar genellikle bunu tercih eder)
    display_position = position + 1 if isinstance(position, int) else position
    display_end_position = end_position + 1 if isinstance(end_position, int) else end_position

    format_args = {
        "position": display_position,
        "end_position": display_end_position,
        "sequence_context": sequence_context or "",
        "length": length,
        "gc_content": gc_content,
        "cpg_oe": cpg_oe, # cpg_oe'yi float olarak alıyoruz
        "enzyme": enzyme,
        "length_aa": length_aa,
        "length_nt": length_nt,
        "frame": frame,
        "strand": strand,
        "algorithm": algorithm,
        **kwargs # Ekstra anahtar kelimeleri yakala
    }

    try:
        # Template içindeki tüm placeholder anahtarlarını al
        placeholders = re.findall(r'\{([^}]+)\}', template)

        # Placeholder'larda format specifier (.2f gibi) varsa onları ayır
        # Örn: "cpg_oe:.2f" -> "cpg_oe"
        placeholder_keys = [key.split(':')[0] for key in placeholders]

        # format_args dictionary'sindeki sadece placeholder'larda bulunan anahtarları içeren yeni bir dictionary oluştur
        # None olanları hariç tutma filtrelemesi kaldırıldı.
        valid_args = {k: format_args[k] for k in placeholder_keys if k in format_args}

        # Şimdi formatlamayı dene
        return template.format(**valid_args)

    except KeyError as e:
        # Eğer gerekli bir anahtar (placeholder'da olan ama format_args'ta olmayan) eksikse buraya düşer
        print(f"Hata: Açıklama formatlamada eksik anahtar: {e}. Template: '{template}', Args: {format_args}")
        return f"{feature_type} (Poz: {display_position}) - Detaylar formatlanamadı (Eksik: {e})."
    except Exception as e:
        # Diğer formatlama veya beklenmedik hatalar için
        print(f"Açıklama formatlama hatası: {e}. Template: '{template}', Args: {format_args}")
        return f"{feature_type} (Poz: {display_position}) - Açıklama oluşturulamadı ({e})."

# Sekans renklendirme ve açıklama için Python fonksiyonu (Frontend JS kullanıyor)
# Bu fonksiyon şu anda doğrudan kullanılmıyor, JS tarafındaki fonksiyonlar daha karmaşık renderlama yapıyor.
# calculate_protein_stats içinde protein HTML dizisini oluştururken kullanılıyor.
def getBaseClass(base, type):
    """Verilen baz/amino asit ve sekans tipine göre CSS sınıfı stringini döndürür."""
    base = base.upper()
    if type == 'DNA':
        if base == 'A': return 'base-A'
        if base == 'T': return 'base-T'
        if base == 'G': return 'base-G'
        if base == 'C': return 'base-C'
        return 'base-Unknown'
    elif type == 'RNA':
        if base == 'A': return 'base-A'
        if base == 'U': return 'base-U'
        if base == 'G': return 'base-G'
        if base == 'C': return 'base-C'
        return 'base-Unknown'
    elif type == 'Protein':
        # Amino asit gruplarına göre sınıf döndür
        if base in ['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'G']: return 'amino-acid-hydrophobic' # Geniş hidrofobik + Glisin/Prolin
        if base in ['S', 'T', 'C', 'Y', 'N', 'Q']: return 'amino-acid-polar'
        if base in ['K', 'R', 'H']: return 'amino-acid-positive'
        if base in ['D', 'E']: return 'amino-acid-negative'
        # Not: Bu Python fonksiyonu sadece Protein HTML çıktısı için sınıf adı döndürüyor.
        # JS tarafındaki renklendirme styles.css değişkenlerini kullanıyor.
        return 'amino-acid-unknown' # X, *, - vb.
    else:
        return 'base-Unknown'


def detect_sequence_type(sequence):
    """Verilen sekansın karakter setine göre türünü (DNA, RNA, Protein, Unknown) belirler."""
    if not sequence:
        return 'Unknown'

    seq_upper = sequence.upper()
    unique_chars = set(seq_upper)

    dna_chars = set("ATGC")
    rna_chars = set("AUGC")
    protein_chars = set("ACDEFGHIKLMNPQRSTVWY*X-") # Standart + Stop + Belirsiz/Boşluk

    is_dna = unique_chars.issubset(dna_chars)
    is_rna = unique_chars.issubset(rna_chars)
    is_protein_candidate = unique_chars.issubset(protein_chars)

    # RNA içeriyorsa (U varsa)
    if 'U' in unique_chars and 'T' not in unique_chars and is_rna:
        return 'RNA'
    # Sadece DNA karakterleri varsa (T olabilir, U olmaz)
    if 'T' in unique_chars and 'U' not in unique_chars and is_dna:
        return 'DNA'
    # Ne T ne U yoksa, ama ATGC'nin alt kümesi ise DNA varsayalım
    if 'T' not in unique_chars and 'U' not in unique_chars and is_dna:
        return 'DNA'
    # Ne T ne U yoksa, ama AUGC'nin alt kümesi ise RNA varsayalım (AUG C gibi kısa diziler)
    # Bu durum yukarıdaki DNA kontrolüyle çakışabilir, öncelik DNA'da
    # if 'T' not in unique_chars and 'U' not in unique_chars and is_rna:
    #     return 'RNA' # Bu durum pek olası değil, ATGC DNA olarak yakalanır

    # DNA veya RNA değilse ve Protein karakter setine uyuyorsa
    if not is_dna and not is_rna and is_protein_candidate:
        return 'Protein'

    # Çok nadir durum: Sadece A, G, C içeriyor olabilir. Bu hem DNA/RNA hem Protein'e uyar.
    # Bu durumda en olası DNA kabul edilir.
    if unique_chars.issubset(set("AGC")):
         return 'DNA'

    return 'Unknown'

def preprocess_sequence(sequence_text):
    """Sekansı temizler (boşlukları kaldırır, büyük harfe çevirir) ve FASTA başlığını atar."""
    if not sequence_text:
        return ""

    lines = sequence_text.strip().splitlines()
    if lines and lines[0].startswith('>'):
        # FASTA başlığını atla, geri kalan satırları birleştir
        sequence = "".join(lines[1:])
    else:
        # Başlık yoksa tüm satırları birleştir
        sequence = "".join(lines)

    # Tüm boşlukları (whitespace) kaldır ve büyük harfe çevir
    cleaned_sequence = re.sub(r'\s+', '', sequence).upper()
    return cleaned_sequence

def find_codons_and_translate_frames(sequence, seq_type):
    """
    DNA/RNA sekansında başlangıç/bitiş kodonlarını bulur ve 6 okuma çerçevesini çevirir.
    Bulunan kodon pozisyonlarını ve çevrilmiş çerçeveleri döndürür.
    Anotasyon listesi döndürür.
    """
    results = {
        "start_codons": [], # 0-based positions
        "stop_codons": [],  # 0-based positions
        "frames": {},       # Translated sequences {'+1': 'M...', ...}
        "annotations": []   # Annotations for codons
    }
    if seq_type not in ['DNA', 'RNA'] or not sequence:
        return results

    seq_obj = Seq(sequence)
    annotations = []

    # Kodonları belirle
    if seq_type == 'DNA':
        start_codon_pattern = "ATG"
        stop_codon_patterns = ["TAA", "TAG", "TGA"]
        rev_comp_seq = seq_obj.reverse_complement()
    else: # RNA
        start_codon_pattern = "AUG"
        stop_codon_patterns = ["UAA", "UAG", "UGA"]
        rev_comp_seq = seq_obj.reverse_complement_rna()

    # Başlangıç ve Bitiş kodonlarını bul (Pozisyonları sakla) ve anotasyon ekle
    for match in re.finditer(start_codon_pattern, sequence):
        pos = match.start()
        results["start_codons"].append(pos)
        annotations.append({
            "type": "start_codon",
            "position": pos,
            "end_position": pos + 2,
            "sequence": match.group(),
            "desc": generate_educational_annotation("start_codon", pos, sequence_context=match.group())
        })

    for stop_pattern in stop_codon_patterns:
        for match in re.finditer(stop_pattern, sequence):
            pos = match.start()
            results["stop_codons"].append(pos)
            annotations.append({
                "type": "stop_codon",
                "position": pos,
                "end_position": pos + 2,
                "sequence": match.group(),
                "desc": generate_educational_annotation("stop_codon", pos, sequence_context=match.group())
            })

    results["start_codons"].sort()
    results["stop_codons"].sort()
    results["annotations"] = annotations # Add codon annotations

    # --- Okuma Çerçevelerini Çevir ---
    sequences_to_process = {
        "+": seq_obj,
        "-": rev_comp_seq
    }

    for strand, current_seq in sequences_to_process.items():
        for frame_num in range(3):
            frame_id = f"{strand}{frame_num + 1}"
            try:
                # Tüm çerçeveyi çevir (stoplar '*' olacak)
                # Biopython translate varsayılan olarak "Standard" tabloyu kullanır (table="Standard")
                # DNA'dan çeviriyorsak (seq_type='DNA'), önce RNA'ya çevirmemiz gerekir.
                if seq_type == 'DNA':
                    temp_rna_seq = current_seq[frame_num:].transcribe()
                    aa_sequence = str(temp_rna_seq.translate(to_stop=False, cds=False))
                else: # seq_type == 'RNA'
                     aa_sequence = str(current_seq[frame_num:].translate(to_stop=False, cds=False))

                # Düzeltme: Burada 'aa_sequence' artık tanımlı

                results["frames"][frame_id] = aa_sequence # Çevrilmiş diziyi frames dictionary'sine ekle

            except Exception as e:
                print(f"Okuma çerçevesi çevirisinde hata ({frame_id}): {e}")
                results["frames"][frame_id] = f"Çeviri hatası: {e}"

    return results

def find_motifs(sequence, seq_type, motifs_dict=None):
    """Verilen sözlükteki regex desenlerine göre sekans içinde motifleri arar."""
    found_motifs = [] # This will store annotation dictionaries
    if not sequence or seq_type not in ['DNA', 'RNA']:
        return found_motifs

    # Genişletilmiş motif desenleri - kullanıcı özel motif sağlamazsa bunlar kullanılır
    # N = [ATGC], R = [AG], Y = [CT], W = [AT], S = [GC], K = [GT], M = [AC], B = [CGT], D = [AGT], H = [ACT], V = [ACG]
    default_motifs = {
        # Temel promotör elementleri (Ökaryotik)
        'TATA-box': r'TATA[AT]A[AT]',          # TATAWAW consensus
        'CAAT-box': r'[GC]CCAAT[CT]',          # GGYCAATCT consensus (Y=C or T)
        'GC-box':   r'GGGCGG|CCGCCC',          # Sp1 binding site consensus (forward/reverse)
        'BRE':      r'[GC][GC][AG]CGCC',       # TFIIB Recognition Element (BRE) upstream of TATA box
        'Inr':      r'[CT][CT]A[ACGT][AT][CT][CT]', # Initiator element (Inr) at transcription start site
        'DPE':      r'[AG]G[AT][CT][GAC][TC]',  # Downstream Promoter Element (DPE)

        # Transkripsiyon düzenleyici elementler
        'enhancer_motif': r'TGA[CG]TCA',       # AP-1 site often found in enhancers
        'silencer_motif': r'CCCTC[ACGT]{5}GAGGG', # Example silencer element pattern
        'AP-1_site':      r'TGA[GC]TCA',       # AP-1 transcription factor binding site
        'NF-kB_site':     r'GGG[AG][AC]TT[CT]CC', # NF-kB binding site consensus (R=A or G, Y=C or T)
        'E-box':          r'CA[ACGT]{2}TG',     # E-box motif (CANNTG)
        'HRE':            r'[AG]CGTG',          # Hipoksi Yanıt Elementi (HRE) - HIF-1 binding site
        'SRE':            r'CC[AT]{6}GG',       # Serum Yanıt Elementi (SRE)
        'CArG-box':       r'CC[AT]{6}GG',       # CArG box (CC(A/T)6GG) - MADS box protein binding site

        # RNA işleme sinyalleri (Ökaryotik)
        'poly_a_signal':  r'AATAAA|ATTAAA',    # Canonical polyadenylation signal and common variant
        'splice_donor':   r'[AC]AGGT[AG]AGT',  # 5' splice site (donor) consensus (MAG|GTRAGT)
        'splice_acceptor':r'[TC]{10,}NCAG[G]?',# 3' splice site (acceptor) consensus (polypyrimidine tract + YAG|G)
        'branch_point':   r'[CT]T[AG]A[CT]',    # Branch point sequence in introns

        # Bakteriyel elementler
        'Shine-Dalgarno': r'AGGAGG[ACGT]{4,8}(ATG|GTG|TTG)', # Shine-Dalgarno sequence (consensus + spacer + start codon)
        'RBS':            r'[AG]{3,9}[ACGT]{4,8}(ATG|GTG|TTG)', # More general Ribosome Binding Site
        'promoter_-35':   r'TTGACA',           # Bacterial -35 promoter element
        'promoter_-10':   r'TATAAT',           # Bacterial -10 promoter element (Pribnow box)
        'rho_independent_terminator': r'GC{6,}[ACGT]{1,10}T{6,}', # Rho-independent terminator (GC-rich stem + T-rich tail)

        # DNA replikasyon elementleri
        'ori_consensus':  r'GATCTNTTNTTTT',    # Example origin consensus (e.g., oriC E.coli fragment)
        'DnaA_box':       r'TT[AT][TC]CCACA',  # DnaA protein binding box consensus (E. coli)
        'autonomously_replicating_sequence': r'[AT]TTTAT[AG]TTT[AT]', # ARS consensus sequence (yeast)

        # Diğer Ökaryotik düzenleyici elementler
        'GATA_motif':     r'[AT]GATA[AG]',     # GATA transcription factor binding site
        'CREB_site':      r'TGACGTCA',         # CREB binding site (cAMP response element - CRE)
        'heat_shock_element': r'[CG][ATCG]GAANNTTC[ATCG][GC]', # Heat shock element (HSE) consensus (N=any)
        'STAT_binding_site': r'TT[CG][CG][GAT]GAA', # STAT transcription factor binding site
        'p53_binding_site': r'[AG][AG][AG]C[AT]{2}G[CT]{2}[CT]', # p53 tumor suppressor binding site
        'NFAT_binding_site': r'GGAAAA',         # NFAT transcription factor binding site
        'ISRE':           r'[AG]GAAA[ACGT]{1,3}GAAA[ACGT]', # İnterferon-Uyarılı Yanıt Elementi

        # Basit Tekrar Örnekleri
        'simple_repeat_AT': r'(AT){5,}',       # Simple AT repeat (at least 5 times)
        'simple_repeat_GC': r'(GC){4,}',       # Simple GC repeat (at least 4 times)
        'simple_repeat_CA': r'(CA){5,}',       # Simple CA repeat (at least 5 times)
        'simple_repeat_GA': r'(GA){5,}',       # Simple GA repeat (at least 5 times)
        'trinucleotide_repeat': r'([ACGT]{3}){4,}', # Trinucleotid Tekrarı (at least 4 times)

        # Telomer ve Sentromer Elementleri
        'telomeric_repeat': r'(TTAGGG){2,}',   # Vertebrate telomeric repeat
        'CENP_B_box': r'TTCGTTGGAAACGGGA',    # Sentromer protein B binding site

        # Metilasyon İlişkili Elementler
        'CpG': r'CG',                         # CpG dinucleotide (potential methylation site)
        'CTCF_binding_site': r'CCGCGNGGNGGCAG', # CTCF insulator protein binding site

        # Mobil Genetik Elementler
        'LTR_consensus': r'TG[AT]{2}[AG]{2}',  # Long Terminal Repeat consensus (simplified)
        'transposon_IR': r'[CT]AGATG[CT]{4}',  # Inverted repeat in transposons (simplified)
    }

    # RNA için U ile değiştir
    if seq_type == 'RNA':
       rna_motifs = {}
       for name, pattern in default_motifs.items():
            # T -> U değişimi ve karakter sınıflarını düzelt
            rna_pattern = pattern
            # Önce karakter sınıflarını düzelt
            rna_pattern = rna_pattern.replace('[AT]', '[AU]')  # W = [AT] -> [AU]
            rna_pattern = rna_pattern.replace('T', 'U')        # Tek başına T'yi U yap
            rna_pattern = rna_pattern.replace('[GT]', '[GU]')  # K = [GT] -> [GU]
            rna_pattern = rna_pattern.replace('[TC]', '[UC]')  # Y = [CT] -> [CU]
            rna_pattern = rna_pattern.replace('[ACT]', '[ACU]') # H = [ACT] -> [ACU]
            rna_pattern = rna_pattern.replace('[AGT]', '[AGU]') # D = [AGT] -> [AGU]
            rna_pattern = rna_pattern.replace('[CGT]', '[CGU]') # B = [CGT] -> [CGU]
            rna_pattern = rna_pattern.replace('[ACGT]', '[ACGU]') # N = [ACGT] -> [ACGU]

            # Start kodon ve stop kodon referanslarını düzelt (eğer regex içinde literal olarak geçiyorsa)
            # Bu regexler zaten U içeriyor, dolayısıyla bu satırlar genellikle gerekmez
            # rna_pattern = rna_pattern.replace('AUG|GUG|UUG', 'AUG|GUG|UUG') # Zaten AUG, GUG, UUG var
            # rna_pattern = rna_pattern.replace('UAA|UAG|UGA', 'UAA|UAG|UGA') # Zaten UAA, UAG, UGA var

            rna_motifs[name] = rna_pattern
       default_motifs = rna_motifs # RNA için güncellenmiş motif setini kullan

    motifs_to_search = motifs_dict if motifs_dict else default_motifs

    for motif_name, pattern in motifs_to_search.items():
        try:
            # re.IGNORECASE ile büyük/küçük harf duyarsız arama
            compiled_pattern = re.compile(pattern, re.IGNORECASE)
            for match in compiled_pattern.finditer(sequence):
                pos = match.start()
                end_pos = match.end() - 1 # Make end position inclusive
                seq_context = match.group()

                # Motif isminin sonuna seq_type eklemeden anotasyon üret
                found_motifs.append({
                    "type": motif_name, # Anotasyon tipi sadece motif adı olacak
                    "position": pos,
                    "end_position": end_pos,
                    "sequence": seq_context,
                    "desc": generate_educational_annotation(motif_name, pos, end_position=end_pos, sequence_context=seq_context)
                })
        except re.error as e:
            print(f"Motif arama hatası - Desen: '{pattern}', Hata: {e}")
        except Exception as e:
            print(f"Beklenmedik motif arama hatası - Desen: '{pattern}', Hata: {e}")


    found_motifs.sort(key=lambda x: x['position'])
    return found_motifs

def calculate_stats(sequence, seq_type):
    """Sekans uzunluğunu, GC içeriğini (DNA/RNA ise) ve nükleotid/amino asit dağılımını hesaplar."""
    stats = {
        "length": 0,
        "gc_content": None,
        "distribution": {},
        "error": None
    }
    if not sequence:
        stats["error"] = "Boş sekans."
        return stats

    try:
        stats["length"] = len(sequence)
        total = stats["length"]

        # Dağılım
        counts = Counter(sequence.upper()) # Büyük harfe çevrilmiş sekans üzerinde sayım yap
        stats["distribution"] = {char: round((count / total) * 100, 2) for char, count in counts.items() if count > 0} # Yüzde olarak, count > 0 olanları al

        # GC İçeriği (Sadece DNA/RNA için)
        if seq_type in ['DNA', 'RNA'] and total > 0:
            try:
                # Biopython's gc_fraction handles non-standard bases by ignoring them
                # It works for both DNA and RNA (counts G+C / total valid bases)
                gc_val = GC_biopython(sequence) # Returns 0.0 to 1.0
                stats["gc_content"] = round(gc_val * 100, 2)
            except Exception as e:
                print(f"GC fraksiyonu hesaplama hatası: {e}")
                stats["gc_content"] = None
                stats["error"] = (stats.get("error", "") + "; " if stats.get("error") else "") + f" GC hesaplama hatası: {e}"
        else:
            stats["gc_content"] = None # Protein veya Unknown için GC içeriği yok

    except Exception as e:
        stats["error"] = f"Genel istatistik hatası: {e}"
        print(f"Hesaplama hatası: {e}\n{traceback.format_exc()}")


    return stats

def find_orfs(sequence, seq_type, min_orf_len_aa=50):
    """
    Verilen DNA/RNA dizisindeki Açık Okuma Çerçevelerini (ORF) bulur.
    Her çerçevede 'M' ile başlayan ve '*' ile biten dizileri arar.
    Anotasyon listesi döndürür.
    """
    orf_annotations = []
    if seq_type not in ['DNA', 'RNA'] or not sequence:
        return orf_annotations

    try:
        seq_obj = Seq(sequence)
        seq_len = len(sequence)

        # Ters tamamlayıcıyı al
        if seq_type == 'DNA':
            rev_comp_seq = seq_obj.reverse_complement()
        else: # RNA
            rev_comp_seq = seq_obj.reverse_complement_rna()

        sequences_to_process = { "+": seq_obj, "-": rev_comp_seq }

        for strand, current_seq_obj in sequences_to_process.items():
            for frame in range(3):
                frame_id = f"{strand}{frame+1}"
                aa_sequence = "" # Her çerçeve için aa_sequence'i tanımla
                try:
                    # Çerçeveyi çevir, stoplar '*' olacak
                    # Biopython translate varsayılan olarak "Standard" tabloyu kullanır (table="Standard")
                    # DNA'dan çeviriyorsak (seq_type='DNA'), önce RNA'ya çevirmemiz gerekir.
                    if seq_type == 'DNA':
                        temp_rna_seq = current_seq_obj[frame:].transcribe()
                        aa_sequence = str(temp_rna_seq.translate(to_stop=False, cds=False))
                    else: # seq_type == 'RNA'
                         aa_sequence = str(current_seq_obj[frame:].translate(to_stop=False, cds=False))

                    # Düzeltme: aa_sequence artık try bloğu içinde tanımlanmış ve döngü öncesine taşınmış

                    # Standard kodon tablosunda sadece ATG Met (M) kodlar. Bakteriyel alternatif başlangıç kodonları
                    # GTG ve TTG de Met kodlayabilir ama standart analizde M'ye odaklanılır.
                    start_codons_aa = {'M'} # Standart Başlangıç AA'sı (Metiyonin)

                    last_start_aa_pos = -1
                    processed_starts = set() # Aynı başlangıçtan birden fazla ORF eklememek için (iç içe M'ler)

                    for i, aa in enumerate(aa_sequence):
                        current_aa_pos = i
                        is_start_aa = (aa in start_codons_aa) # Belirlenen başlangıç AA'larından biri mi?

                        # ORF başlangıcını bul
                        if is_start_aa:
                            last_start_aa_pos = current_aa_pos
                            # Bu başlangıç AA'sından sonraki ilk stop'u bul
                            first_stop_after_start_aa_pos = -1
                            for j in range(last_start_aa_pos + 1, len(aa_sequence)):
                                if aa_sequence[j] == '*':
                                    first_stop_after_start_aa_pos = j
                                    break # İlk stop bulundu

                            # Eğer bir stop kodonu bulunduysa
                            if first_stop_after_start_aa_pos != -1:
                                # ORF'un AA dizisi (başlangıç dahil, stop hariç)
                                orf_aa_seq = aa_sequence[last_start_aa_pos : first_stop_after_start_aa_pos]
                                orf_len_aa = len(orf_aa_seq)

                                # Minimum uzunluk kontrolü
                                if orf_len_aa >= min_orf_len_aa:
                                    # Nükleotid pozisyonlarını hesapla (0-tabanlı)
                                    # AA pozisyonundan NT pozisyonuna: nt_pos = frame_offset + aa_pos * 3
                                    # Başlangıç kodonunun başlangıç pozisyonu
                                    orf_start_nt_0based = frame + last_start_aa_pos * 3
                                    # Stop kodonunun başlangıç pozisyonu
                                    stop_codon_start_nt_0based = frame + first_stop_after_start_aa_pos * 3
                                    # Stop kodonunun son bazı (dahil)
                                    orf_end_nt_0based = stop_codon_start_nt_0based + 2

                                    # Nükleotid dizisini al (start dahil, stop dahil)
                                    # Bu slice'ı orijinal (işlenmemiş) sekans üzerinden yapmalıyız
                                    # Ters strand ise pozisyonları orijinal sekansa göre çevirmemiz gerekir
                                    actual_nt_start = orf_start_nt_0based
                                    actual_nt_end = orf_end_nt_0based

                                    if strand == "-":
                                        # Ters strandde hesaplanan 0-tabanlı pozisyonları orijinal sekansın 0-tabanlı pozisyonlarına çevir
                                        # current_seq_obj ters tamamlayıcı diziyi temsil ediyordu
                                        # current_seq_obj'deki index `k` -> orijinal seq objesindeki index `seq_len - 1 - k`
                                        # orf_start_nt_0based (current_seq_obj üzerindeki start) -> orig_seq_obj üzerindeki END
                                        # orf_end_nt_0based (current_seq_obj üzerindeki end)   -> orig_seq_obj üzerindeki START
                                        orig_start_nt = seq_len - 1 - orf_end_nt_0based
                                        orig_end_nt   = seq_len - 1 - orf_start_nt_0based

                                        actual_nt_start = orig_start_nt
                                        actual_nt_end = orig_end_nt
                                        # Sağdan sola okuduğumuz için aslında başlangıç > bitiş gibi görünebilir.
                                        # Anotasyonlarda start < end olması görselleştirme için daha iyi olabilir.
                                        # Orijinal sekanstaki en küçük pozisyonu start, en büyük pozisyonu end yapalım.
                                        orf_start_nt_for_anno = min(orig_start_nt, orig_end_nt)
                                        orf_end_nt_for_anno = max(orig_start_nt, orig_end_nt)


                                    else: # Düz strand
                                         orf_start_nt_for_anno = actual_nt_start
                                         orf_end_nt_for_anno = actual_nt_end


                                    orf_len_nt = (orf_end_nt_for_anno - orf_start_nt_for_anno + 1) # NT uzunluğu pozisyon farkı + 1

                                    orf_nt_seq = ""
                                    try:
                                         # Nükleotid dizisini almak için orijinal seq objesini kullan
                                         # slice [start : end + 1] (Python slicingde end dahil değil, o yüzden +1)
                                         orf_nt_seq = str(seq_obj[actual_nt_start : actual_nt_end + 1])

                                    except IndexError:
                                        orf_nt_seq = "[Hata: Dizi Sınırları Aşıldı]"
                                    except Exception as slice_err:
                                         print(f"ORF NT slice error: {slice_err}")
                                         orf_nt_seq = "[Hata: Slice Başarısız]"

                                    # Eğer bu başlangıç pozisyonundan (orijinal dizi üzerindeki) zaten bir ORF eklediysek, atla
                                    # (Birden çok ORF aynı ATG'de başlayabilir ama sadece ilkini listelemek isteyebiliriz)
                                    if orf_start_nt_for_anno in processed_starts:
                                        continue
                                    processed_starts.add(orf_start_nt_for_anno)


                                    orf_annotations.append({
                                        "type": "ORF",
                                        "position": orf_start_nt_for_anno, # Anotasyon için düzenlenmiş başlangıç (küçük olan)
                                        "end_position": orf_end_nt_for_anno, # Anotasyon için düzenlenmiş bitiş (büyük olan)
                                        "frame": frame_id,
                                        "strand": strand,
                                        "length_aa": orf_len_aa,
                                        "length_nt": orf_len_nt,
                                        "protein_sequence": orf_aa_seq,
                                        "nucleotide_sequence": orf_nt_seq, # Start ve Stop kodonunu içeren nt dizisi (orijinal yönde)
                                        "desc": generate_educational_annotation(
                                            "ORF", orf_start_nt_for_anno, end_position=orf_end_nt_for_anno,
                                            length_aa=orf_len_aa, length_nt=orf_len_nt, frame=frame_id, strand=strand
                                        )
                                    })
                                # Başarılı bir ORF bulduktan sonra, bu M için aramayı bitir.
                                # Overlapping M'lerden başlayanları listelemeye devam etmek istiyorsak bu satırı YORUM SATIRI yapın:
                                # last_start_aa_pos = -1 # Resetlemeye gerek yok, döngü devam edecek (Zaten for döngüsü devam ediyor)

                except Exception as e:
                    # Düzeltme: Eğer çeviri sırasında hata oluşursa (temp_rna_seq veya translate)
                    # aa_sequence tanımlanmayabilir. Bu durumda döngüye girmeden hatayı işle.
                    print(f"ORF bulma hatası - Çerçeve {frame_id}: Çeviri sırasında hata veya aa_sequence erişimi: {e}\n{traceback.format_exc()}")
                    # Bu çerçeve için ORF bulunamadı olarak devam eder.


    except Exception as e:
        print(f"Genel ORF bulma hatası: {e}\n{traceback.format_exc()}")
        # Hata durumunda boş liste döndürülmüş olacak

    # ORF'ları uzunluğa göre büyükten küçüğe sırala
    orf_annotations.sort(key=lambda x: x.get('length_aa', 0), reverse=True)
    return orf_annotations

def find_restriction_sites(sequence, seq_type, enzyme_list=None):
    """Verilen sekans üzerinde restriksiyon enzimi kesim alanlarını bulur ve anotasyon listesi döndürür."""
    site_annotations = []
    if seq_type != 'DNA' or not sequence: # Genellikle DNA'da yapılır
        return site_annotations

    # CommOnly: Yaygın kullanılan enzimler seti
    # Veya kullanıcıdan bir liste alınabilir.
    # Eğer enzyme_list None ise varsayılanı kullan
    if enzyme_list is None:
        # Temel, yaygın birkaç enzimle başlayalım
        # Biopython CommOnly RestrictionBatch'i de kullanılabilir, bu daha kapsamlıdır.
        # batch = RestrictionBatch(CommOnly)
        # Veya belirli bir liste:
        enzyme_list_to_use = ['EcoRI', 'BamHI', 'HindIII', 'NotI', 'XbaI', 'SpeI', 'PstI', 'SacI', 'KpnI']
    else:
         # Kullanıcıdan gelen listeyi kullan, boş veya geçersiz olabilir
         # Geçersiz isimler Biopython tarafından göz ardı edilebilir, bu nedenle batch boş kalabilir
         enzyme_list_to_use = enzyme_list

    try:
        seq_obj = Seq(sequence)
        # Seçili enzim listesi ile RestrictionBatch oluştur
        # Biopython'da enzim adları büyük/küçük harf duyarlı olabilir, emin olmak için listeyi normalize edebiliriz
        # Ancak RestrictionBatch genellikle popüler isimleri tanır.
        batch = RestrictionBatch(enzyme_list_to_use) # Kullanıcının verdiği veya varsayılan listeyi kullan

        # Eğer batch boşsa (geçersiz enzim adı verildi veya liste boş), analiz yapma
        if not batch:
            # Hata mesajı döndür ve çık
            msg = f"Belirtilen restriksiyon enzimleri ({enzyme_list_to_use}) bulunamadı veya geçersiz."
            print(f"Restriksiyon analizi hatası: {msg}")
            return [{"type": "error", "message": msg}]


        analysis = Analysis(batch, seq_obj, linear=True) # Dizinin linear olduğunu varsayalım

        # .search() metodunu kullanarak site objelerine ve pozisyonlara ulaşabiliriz
        # search() 1-tabanlı pozisyonlar döndürür
        results_with_sites = analysis.search(seq_obj) # {EnzymeObj: [pos1based_1, pos1based_2, ...]}


        for enzyme_obj, positions1based in results_with_sites.items():
             enzyme_name = str(enzyme_obj)
             # enzyme_obj'nin site özniteliği var mı kontrol et
             recognition_site = enzyme_obj.site if hasattr(enzyme_obj, 'site') else 'N/A'
             site_length = len(recognition_site) if recognition_site != 'N/A' else 0

             for pos1based in positions1based: # pos1based burada sadece pozisyon (int)
                 pos0based = pos1based - 1
                 # Tanıma bölgesinin sonu (dahil, 0-tabanlı)
                 # Site uzunluğu kadar ileri git, ama son index dahil olduğu için -1
                 end_pos0based = pos0based + site_length - 1 if site_length > 0 else pos0based

                 # Anotasyon pozisyonlarının dizinin sınırları içinde olduğundan emin ol
                 if pos0based >= 0 and pos0based < len(sequence): # Başlangıç pozisyonu geçerli mi?
                      # Bitiş pozisyonu başlangıçtan küçük olabilir mi? Hayır, recognition site uzunluğu > 0 ise end >= start.
                      # Sadece end_pos0based'in dizi içinde olup olmadığını kontrol edelim
                      if end_pos0based >= 0 and end_pos0based < len(sequence):
                          site_annotations.append({
                              "type": "restriction_site",
                              "position": pos0based,
                              "end_position": end_pos0based, # Tanıma bölgesinin sonu (dahil)
                              "enzyme": enzyme_name,
                              "sequence": recognition_site, # Tanıma dizisi
                              "desc": generate_educational_annotation("restriction_site", pos0based, enzyme=enzyme_name, sequence_context=recognition_site)
                          })
                      else:
                           # Bitiş dizi sınırları dışında ama başlangıç içinde. Sadece başlangıcı işaretle.
                           # Veya bu durumu farklı bir anotasyon tipiyle göster. Şimdilik sadece başlangıcı işaretleyelim.
                           # Anotasyon tanımımız aralık bekliyor, bu yüzden end_position'ı başlangıçla aynı yapalım.
                           site_annotations.append({
                              "type": "restriction_site",
                              "position": pos0based,
                              "end_position": pos0based, # Sadece başlangıcı işaretle
                              "enzyme": enzyme_name + " (Kısmi)", # Kısmi olduğunu belirt
                              "sequence": recognition_site, # Tanıma dizisi
                              "desc": generate_educational_annotation("restriction_site", pos0based, enzyme=enzyme_name + " (Kısmi)", sequence_context=recognition_site) + " (Alan dizi sınırları dışında bitiyor)"
                           })
                           print(f"Uyarı: Restriksiyon alanı {enzyme_name} @ {pos1based} ({pos0based}) dizi sınırları dışında bitiyor. Başlangıç işaretlendi.")

                 else:
                     print(f"Uyarı: Restriksiyon alanı {enzyme_name} @ {pos1based} ({pos0based}) dizi sınırları dışında. Atlanıyor.")


    except ValueError as ve: # Örneğin listede olmayan enzim adı veya Biopython'ın ValueError'u
         msg = f"Restriksiyon analizi hatası (Geçersiz enzim adı veya Biopython problemi): {ve}"
         print(msg)
         # Hata durumunda hata objesi içeren liste döndür (app.py bunu yakalar)
         return [{"type": "error", "message": msg}]
    except Exception as e:
        msg = f"Restriksiyon analizi sırasında beklenmedik hata: {e}"
        print(msg + f"\n{traceback.format_exc()}")
        return [{"type": "error", "message": msg}]

    site_annotations.sort(key=lambda x: x['position'])
    return site_annotations

def calculate_cpg_oe(sequence_window):
    """Verilen penceredeki CpG gözlenen/beklenen oranını hesaplar."""
    seq_len = len(sequence_window)
    if seq_len == 0:
        return 0.0

    # Sadece ATGC sayalım, diğer karakterleri yok sayalım
    valid_chars = "ATGC"
    cleaned_window = ''.join(c for c in sequence_window.upper() if c in valid_chars)
    cleaned_len = len(cleaned_window)

    if cleaned_len < 2: # En az 2 baz lazım
        return 0.0

    c_count = cleaned_window.count('C')
    g_count = cleaned_window.count('G')
    cpg_count = 0

    # CpG dinükleotidlerini doğru şekilde say (örtüşen CG'leri saymak için)
    for j in range(cleaned_len - 1):
        if cleaned_window[j:j+2] == 'CG':
            cpg_count += 1

    # Bekleneni hesapla
    # Eğer C veya G yoksa veya toplam uzunluk yetersizse O/E tanımsızdır (veya 0 kabul edilebilir)
    if c_count == 0 or g_count == 0:
        return 0.0

    # Formül: O/E = (Observed CpG * Total_len) / (Count C * Count G)
    # Bu formül Gardiner-Garden & Frommer'ın kullandığı formülle uyumludur.
    expected_cpg = (c_count * g_count) / cleaned_len # Beklenen CG sayısı

    if expected_cpg == 0:
         return 0.0 # Bölme hatasını önle

    cpg_oe = cpg_count / expected_cpg # Gözlenen / Beklenen oranı

    return cpg_oe

def find_cpg_islands(sequence, seq_type, algorithm='gardiner', min_length=None, window_size=None, step_size=1, min_gc_override=None, min_oe_override=None):
    """
    CpG adalarını tespit eder. Farklı algoritmalar kullanılabilir. Anotasyon listesi döndürür.

    Parametreler:
        sequence (str): DNA dizisi
        seq_type (str): Dizi tipi ('DNA' olmalı)
        algorithm (str): 'gardiner', 'takai', 'ucsc', 'irizarry' (farklı parametre setleri)
        min_length (int): Ada için minimum uzunluk (bp) - Varsayılan algoritmaya göre değişir, override edilebilir.
        window_size (int): Kaydırma penceresi boyutu (bp) - Varsayılan algoritmaya göre değişir, override edilebilir.
        step_size (int): Pencere kaydırma adımı (bp) - 1 genellikle en hassas
        min_gc_override (float): Kullanıcının belirlediği min GC% (algoritma varsayılanını geçersiz kılar)
        min_oe_override (float): Kullanıcının belirlediği min O/E (algoritma varsayılanını geçersiz kılar)

    Algoritma Kriterleri ve Varsayılanları:
    - Gardiner-Garden & Frommer (1987): min_gc=50%, min_oe=0.6, min_length=200 bp, window_size=100 bp (Yaygın kullanımda 100 veya 200)
    - Takai & Jones (2002): min_gc=55%, min_oe=0.65, min_length=500 bp, window_size=200 bp
    - UCSC (Approx): min_gc=50%, min_oe=0.6, min_length=200 bp (Genellikle HMM kullanır, burada basit eşikleme)
    - Irizarry et al. (2009): min_gc=50%, min_oe=0.6, min_length=200 bp, window_size=100 bp
    """
    island_annotations = []
    if seq_type != 'DNA' or not sequence:
        return island_annotations

    seq_len = len(sequence)
    sequence_upper = sequence.upper() # Emin olmak için

    # Algoritma varsayılanlarını ayarla (min_length ve window_size artık buradan gelecek varsayılanları)
    # Eğer frontend'den gelmezse None olacak, o zaman buradaki varsayılanları kullanacağız.
    default_params = {
        'gardiner': {'min_gc': 50.0, 'min_oe': 0.6, 'min_length': 200, 'window_size': 100}, # Daha sık kullanılan window size 100
        'takai':    {'min_gc': 55.0, 'min_oe': 0.65, 'min_length': 500, 'window_size': 200},
        'ucsc':     {'min_gc': 50.0, 'min_oe': 0.6, 'min_length': 200, 'window_size': 200}, # UCSC aslında 200 bp min length kullanıyor gibi
        'irizarry': {'min_gc': 50.0, 'min_oe': 0.6, 'min_length': 200, 'window_size': 100}
    }

    if algorithm not in default_params:
        algorithm = 'gardiner' # Geçersizse varsayılana dön

    params = default_params[algorithm]

    # Kullanıcı override'larını veya algoritma varsayılanlarını uygula
    min_gc_thresh = min_gc_override if min_gc_override is not None else params['min_gc']
    min_oe_thresh = min_oe_override if min_oe_override is not None else params['min_oe']
    current_min_length = min_length if min_length is not None else params['min_length']
    current_window_size = window_size if window_size is not None else params['window_size']

    # Minimum uzunluk veya pencere boyutu sekans uzunluğundan büyükse çık
    if seq_len < current_window_size or seq_len < current_min_length:
        return island_annotations

    # 1. Adım: Kriterleri sağlayan pencereleri bul
    candidate_indices = []
    # Pencere başlangıcı `i` için döngü. Son pencere `seq_len - current_window_size` konumunda başlar.
    for i in range(0, seq_len - current_window_size + 1, step_size):
        window = sequence_upper[i : i + current_window_size]

        # GC içeriğini hesapla (sadece ATGC kullanarak)
        gc_content = GC_biopython(window) * 100

        # Gözlenen/Beklenen CpG oranını hesapla
        cpg_oe_val = calculate_cpg_oe(window) # Fonksiyon artık yuvarlama yapmıyor

        # Şartları kontrol et
        if gc_content >= min_gc_thresh and cpg_oe_val >= min_oe_thresh:
            candidate_indices.append(i) # Bu pencerenin başlangıç indeksini ekle

    if not candidate_indices:
        return island_annotations # Hiç aday pencere yoksa

    # 2. Adım: Bitişik veya çok yakın aday pencereleri birleştirerek adaları oluştur
    # İki aday pencere arasındaki maksimum boşluk (gap) ne kadar olabilir?
    # Genellikle pencere boyutunun yarısı veya sabit küçük bir değer (örn. 20-50 bp)
    # Kullanılan algoritmaya göre değişebilir. Bazı tanımlar step_size'a eşit gap'e izin verir.
    # Basit birleştirme: Eğer bir sonraki aday penceresi, mevcut pencerenin bitişinden
    # belirlenen bir "izin verilen boşluk" (örneğin 20 bp gibi küçük bir değer) kadar uzaktaysa birleştir.
    # Yaygın birleştirme kuralı: Eğer pencereler arasındaki boşluk <= 20 bp ise birleştir.
    merge_gap_threshold = 20 # Yaygın kabul gören maksimum boşluk (bp)

    merged_islands_raw = []
    if not candidate_indices: return island_annotations

    current_island_start = candidate_indices[0]
    # current_island_end, birleştirilmiş adanın gerçek son pozisyonunu tutacak
    current_island_end = candidate_indices[0] + current_window_size - 1

    for i in range(1, len(candidate_indices)):
        prev_window_end = candidate_indices[i-1] + current_window_size - 1
        current_window_start = candidate_indices[i]

        # Eğer mevcut pencerenin başlangıcı, bir önceki pencerenin bitişinden hemen sonra
        # (aradaki boşluk <= merge_gap_threshold ise) adayı uzat.
        # Boşluk = current_window_start - (prev_window_end + 1)
        gap = current_window_start - (prev_window_end + 1)

        if gap <= merge_gap_threshold:
            # Birleştirilmiş adanın bitişini, mevcut pencerenin bitişi ile güncelle
            current_island_end = current_window_start + current_window_size - 1
        else:
            # Ada bitti, kaydet ve yeni bir ada başlat
            # Adanın bitişi son pencerenin bitişiydi
            merged_islands_raw.append({'start': current_island_start, 'end': current_island_end})
            # Yeni ada başlangıcı ve bitişi
            current_island_start = current_window_start
            current_island_end = current_window_start + current_window_size - 1

    # Döngüden sonra son adayı da ekle
    merged_islands_raw.append({'start': current_island_start, 'end': current_island_end})


    # 3. Adım: Birleştirilmiş adaların uzunluğunu kontrol et ve detayları hesapla/ekle
    for island in merged_islands_raw:
        start = island['start']
        end = island['end']
        length = end - start + 1

        # Adanın minimum uzunluk kriterini sağlıyor mu kontrol et
        if length >= current_min_length:
            # Adanın gerçek dizisini al
            island_seq = sequence_upper[start : end + 1]
            # Ada için GC ve O/E'yi yeniden hesapla (bu kez tüm ada için)
            final_gc = round(GC_biopython(island_seq) * 100, 2)
            final_oe = calculate_cpg_oe(island_seq) # Yuvarlama yapma, sonradan yaparız

            # Ada hala genel GC ve O/E kriterlerini sağlıyor mu diye kontrol et
            # Gardiner/UCSC için bu kontrol genellikle istenir. Takai daha sıkı.
            # Genel bir kontrol olarak ekleyelim:
            if final_gc >= min_gc_thresh and final_oe >= min_oe_thresh:

                island_annotations.append({
                    "type": "cpg_island",
                    "position": start,
                    "end_position": end,
                    "length": length,
                    "gc_content": final_gc, # Zaten yuvarlandı
                    "cpg_oe": round(final_oe, 3), # O/E'yi burada yuvarlayalım
                    "algorithm": algorithm, # Kullanılan algoritma adını sakla
                    "desc": generate_educational_annotation(
                        "cpg_island", start, end_position=end, length=length,
                        gc_content=final_gc, cpg_oe=final_oe, algorithm=algorithm
                    )
                })
            # else:
                # print(f"Ada atlandı: {start}-{end}, Uz: {length}, GC: {final_gc}%, O/E: {final_oe:.2f}")


    island_annotations.sort(key=lambda x: x['position'])
    return island_annotations


def calculate_protein_stats(protein_sequence):
    """
    Verilen protein dizisi için temel istatistikleri hesaplar ve okunabilir bir formatta sunar.
    Ayrıca eğitimsel fonksiyon ipuçlarını da döndürür.
    """
    stats = { "error": None }
    if not protein_sequence or not isinstance(protein_sequence, str):
        stats["error"] = "Geçerli bir protein dizisi değil."
        stats["potential_function_hints"] = ["Protein dizisi sağlanamadı veya geçersiz."]
        return stats

    # Geçersiz karakterleri temizle (ProtParam bunları kabul etmez)
    # Standart 20 AA dışındakileri temizleyelim. X, *, - kabul edilmez.
    valid_aa_chars = "ACDEFGHIKLMNPQRSTVWY"
    protein_sequence_cleaned = ''.join(aa for aa in protein_sequence.upper() if aa in valid_aa_chars)

    if not protein_sequence_cleaned:
        stats["error"] = "Analiz edilecek geçerli amino asit bulunamadı (Temizleme sonrası boş)."
        stats["potential_function_hints"] = ["Temizleme sonrası protein dizisi boş."]
        return stats

    try:
        # ProtParam objesi oluştur
        # Biopython >= 1.80 Seq objesi bekler
        try:
             analyzer = ProtParam.ProteinAnalysis(Seq(protein_sequence_cleaned))
        except TypeError:
             # Eski Biopython sürümleri string kabul edebilir
             analyzer = ProtParam.ProteinAnalysis(protein_sequence_cleaned)


        # Amino asit grupları (fizikokimyasal özelliklerine göre)
        aa_groups = {
            "hydrophobic": ("A", "V", "I", "L", "M", "F", "Y", "W"), # Hidrofobik (Nonpolar)
            "polar": ("S", "T", "N", "Q"),          # Polar (Nötr)
            "positive": ("R", "H", "K"),            # Pozitif yüklü (Bazik)
            "negative": ("D", "E"),                 # Negatif yüklü (Asidik)
            "special": ("C", "G", "P")              # Özel (Sistein, Glisin, Prolin)
        }

        # Amino asit kısaltmaları ve tam adları
        aa_names = {
            'A': 'Alanin', 'C': 'Sistein', 'D': 'Aspartik Asit', 'E': 'Glutamik Asit',
            'F': 'Fenilalanin', 'G': 'Glisin', 'H': 'Histidin', 'I': 'İzolösin',
            'K': 'Lizin', 'L': 'Lösin', 'M': 'Metiyonin', 'N': 'Asparajin',
            'P': 'Prolin', 'Q': 'Glutamin', 'R': 'Arjinin', 'S': 'Serin',
            'T': 'Treonin', 'V': 'Valin', 'W': 'Triptofan', 'Y': 'Tirozin'
        }

        # Amino asit renk kodları (fizikokimyasal özelliklerine göre) - Frontend için
        # getBaseClass fonksiyonu bu renk sınıflarını döndürüyor
        # Bu dictionary sadece bilgi amaçlı, JS'de kullanılacak ve doğrudan döndürülmüyor
        # aa_colors = { # Bu dictionary sadece bilgi amaçlı, JS'de kullanılacak
        #     'A': '#3498db', 'V': '#2980b9', 'L': '#1f618d', 'I': '#1a5276',
        #     'M': '#154360', 'F': '#0e2f44', 'W': '#0a2229', 'Y': '#9b59b6',
        #     'S': '#2ecc71', 'T': '#27ae60', 'N': '#229954', 'Q': '#1e8449',
        #     'K': '#e74c3c', 'R': '#c0392b', 'H': '#a93226',
        #     'D': '#f39c12', 'E': '#d35400',
        #     'C': '#95a5a6', 'G': '#7f8c8d', 'P': '#34495e'
        # }


        # Protein dizisini formatla (10'ar aa gruplar halinde, pozisyon numaraları ile) - Frontend için
        formatted_sequence = ""
        line_length_aa = 10 # Her grupta kaç AA
        group_spacing = "  " # Gruplar arası boşluk
        line_spacing = "\n" # Satırlar arası boşluk
        positions_padding = 6 # Pozisyon numarası için dolgu

        lines = []
        # Her 6 AA grubunda (toplam 60 AA) bir satır sonu
        for i in range(0, len(protein_sequence_cleaned), line_length_aa * 6):
            line_chunk = protein_sequence_cleaned[i : i + line_length_aa * 6]
            formatted_line = f"{str(i+1).rjust(positions_padding)}  " # Pozisyon numarası

            for j in range(0, len(line_chunk), line_length_aa):
                group_seq = line_chunk[j : j + line_length_aa]
                # İçteki 5'lik gruplar arası boşluk (Biopython ProtParam çıktısına benzer format)
                sub_groups = [group_seq[k:k+5] for k in range(0, len(group_seq), 5)]
                formatted_line += " ".join(sub_groups) + group_spacing

            lines.append(formatted_line.strip()) # Satır sonundaki fazla boşluğu kaldır

        formatted_sequence = line_spacing.join(lines)


        # HTML formatında renkli dizi gösterimi - Frontend için (Bu kodun Python'da olması yerine JS'de olması daha mantıklı)
        # Ancak şimdilik Python'daki getBaseClass'ı kullanarak HTML snippet'i oluşturalım
        html_sequence = "<div class='protein-sequence-html'>"
        html_line_length = 60 # Toplam karakter (pozisyon dahil değil)
        html_group_length = 10 # Her grupta 10 AA

        for i in range(0, len(protein_sequence_cleaned), html_line_length):
            line_seq = protein_sequence_cleaned[i:i+html_line_length]
            html_sequence += f"<div class='sequence-line'><span class='position-number'>{i+1}</span> "
            for j in range(0, len(line_seq), html_group_length):
                group_seq = line_seq[j:j+html_group_length]
                html_sequence += "<span class='aa-group'>"
                for aa in group_seq:
                    # getBaseClass fonksiyonunu burada Python'da çağırıp sınıf adını alıyoruz
                    css_class = getBaseClass(aa, 'Protein')
                    title_text = aa_names.get(aa, aa) # Tam ad title'da gösterilebilir
                    html_sequence += f"<span class='aa {css_class}' title='{title_text}'>{aa}</span>"
                html_sequence += "</span> " # Grup sonu boşluk
            html_sequence += f" <span class='position-number'>{i + len(line_seq)}</span></div>" # Satır sonu pozisyon
        html_sequence += "</div>"


        # Amino asit gruplarına göre dağılımı hesapla
        aa_group_counts = {group: 0 for group in aa_groups}
        for aa in protein_sequence_cleaned:
            for group, aas in aa_groups.items():
                if aa in aas:
                    aa_group_counts[group] += 1
                    break # Bir gruba aitse diğerlerine bakma

        # Grup yüzdelerini hesapla
        total_aa = len(protein_sequence_cleaned)
        aa_group_percent = {group: round((count / total_aa) * 100, 2)
                           for group, count in aa_group_counts.items()} if total_aa > 0 else {}

        # Protein özellikleri sınıflandırması
        instability_index = analyzer.instability_index()
        stability = "Kararlı (Stabil)" if instability_index < 40 else "Kararsız (Instabil)"

        gravy_value = analyzer.gravy()
        hydropathy = "Hidrofilik (Suda Çözünür Eğilimli)" if gravy_value < 0 else "Hidrofobik (Yağda Çözünür/Membran Eğilimli)"

        # Amino asit yüzdelerini al ve isimlendir
        aa_percent_raw = analyzer.get_amino_acids_percent()
        aa_percent_named = {aa_names.get(aa, aa): round(perc * 100, 2) for aa, perc in aa_percent_raw.items() if perc > 0} # Sadece >0 olanları al

        # Amino asit sayılarını da ekle
        aa_counts = {aa: protein_sequence_cleaned.count(aa) for aa in valid_aa_chars if protein_sequence_cleaned.count(aa) > 0}
        aa_counts_named = {aa_names.get(aa, aa): count for aa, count in aa_counts.items()}

        # Protein özellikleri açıklamaları - Frontend için
        property_explanations = {
            "molecular_weight": "Proteinin dalton (Da) cinsinden moleküler ağırlığı. 1 Da yaklaşık olarak bir hidrojen atomunun kütlesine eşittir.",
            "isoelectric_point": "Proteinin net yükünün sıfır olduğu pH değeri. Bu pH'da protein çözünürlüğü en düşüktür.",
            "instability_index": "Proteinin in vitro kararlılığının bir ölçüsü. 40'tan düşük değerler kararlı proteinleri gösterir.",
            "aromaticity": "Proteindeki aromatik amino asitlerin (Phe, Tyr, Trp) oranı.",
            "gravy": "Grand Average of Hydropathy. Pozitif değerler hidrofobik, negatif değerler hidrofilik proteinleri gösterir."
        }

        # Eğitimsel fonksiyon ipuçlarını burada hesapla
        potential_function_hints = get_protein_function_hints(protein_sequence_cleaned, analyzer)


        stats = {
            "length": analyzer.length,
            "molecular_weight": round(analyzer.molecular_weight(), 2),
            "isoelectric_point": round(analyzer.isoelectric_point(), 2),
            "amino_acid_percent_named": aa_percent_named, # İsimlendirilmiş yüzdeler
            "amino_acid_counts": aa_counts_named, # Amino asit sayıları
            "amino_acid_groups": aa_group_percent,  # Amino asit grup dağılımları
            "aromaticity": round(analyzer.aromaticity(), 3),
            "instability_index": round(instability_index, 2),
            "stability_class": stability,
            "gravy": round(gravy_value, 3), # Hidrofobisite
            "hydropathy_class": hydropathy,
            "secondary_structure_fraction": {
                 struct: round(perc * 100, 1) for struct, perc in zip(['Helix', 'Turn', 'Sheet'], analyzer.secondary_structure_fraction())
            },
            "formatted_sequence": formatted_sequence, # Kullanıcı arayüzü için formatlanmış dizi (pre tag içine)
            "html_sequence": html_sequence, # HTML formatında renkli dizi (span'lerle, Python getBaseClass kullanarak)
            "full_sequence": protein_sequence_cleaned, # Tam, temizlenmiş dizi
            "property_explanations": property_explanations, # Özellik açıklamaları (Frontend için)
            # "aa_colors": aa_colors, # Amino asit renk kodları (Frontend için - artık JS'de styles.css değişkenlerinden alınıyor)
            "potential_function_hints": potential_function_hints # Eğitimsel ipuçları
        }
        return stats
    except Exception as e:
        print(f"Protein analizi hatası: {e}\n{traceback.format_exc()}")
        stats["error"] = f"Protein analizi sırasında hata: {e}"
        stats["potential_function_hints"] = [f"Protein analizi sırasında ipuçları alınamadı: {e}"]
        return stats

def get_protein_function_hints(protein_sequence, analyzer):
    """
    (EĞITIMSEL) Protein dizisinin özelliklerine bakarak olası fonksiyon ipuçları verir.
    NOT: Bu fonksiyon sadece eğitimsel amaçlıdır ve gerçek protein fonksiyon tahmini için ASLA kullanılmamalıdır.
    Bu fonksiyon calculate_protein_stats tarafından çağrılacaktır.
    """
    hints = []
    # Analiz için yeterli protein dizisi ve analyzer objesi olup olmadığını kontrol et
    if not protein_sequence or not analyzer:
        hints.append("Analiz için yeterli protein verisi veya araç objesi yok.")
        return hints
     # Minimum bir uzunluk kontrolü ekleyelim
    if len(protein_sequence) < 10:
         hints.append("Çok kısa protein dizisi, anlamlı ipuçları verilemiyor.")
         return hints


    # --- Basit Motif Kontroleri ---
    # Cistein çiftleri (metal bağlama, disülfit)
    if re.search(r"C..C", protein_sequence):
        hints.append("Cys-X-X-Cys (C..C) motifi: Disülfit bağı oluşumu veya metal iyonu (örn: çinko) bağlama potansiyeli.")
    if re.search(r"CXXC", protein_sequence):
        hints.append("Cys-X-X-Cys (CXXC) motifi: Redoks aktif bölge potansiyeli (örn: tioredoksin benzeri proteinler).")
    # Yüksek Prolin içeriği
    proline_count = protein_sequence.count('P')
    if len(protein_sequence) > 0 and proline_count / len(protein_sequence) > 0.1:
         hints.append(f"Yüksek Prolin içeriği (%{(proline_count/len(protein_sequence))*100:.1f}): Yapısal olarak sert veya bükülgen olmayan bölgeler, protein-protein etkileşimleri (örn: SH3 domainleri ile etkileşim).")
    # Yüksek Glisin içeriği
    glycine_count = protein_sequence.count('G')
    if len(protein_sequence) > 0 and glycine_count / len(protein_sequence) > 0.1:
         hints.append(f"Yüksek Glisin içeriği (%{(glycine_count/len(protein_sequence))*100:.1f}): Esnek bölgeler, sıkı dönüşler veya yapısal olmayan bölgeler (Intrinsically Disordered Regions).")
    # ER Retansiyon Sinyalleri
    if protein_sequence.endswith("KDEL") or protein_sequence.endswith("HDEL"):
        hints.append("KDEL/HDEL C-terminal motifi: Endoplazmik Retikulum'da (ER) kalma sinyali olabilir.")
    # Nükleer Lokalizasyon Sinyali (Çok Basit)
    # Klasik NLS: K(K/R)X(K/R) veya PXXKKR benzeri
    if re.search(r"[KR]{2,}[^P]", protein_sequence) or re.search(r"P[KR]{2,}[^P]", protein_sequence) or re.search(r"[KR]{4,}", protein_sequence):
        hints.append("Bazik amino asit kümesi veya PXXKR motifi benzeri: Nükleer lokalizasyon sinyali (NLS) potansiyeli.")
    # Hücre İçi Hedeflenme/Sıralama Sinyalleri (Çok Basit)
    if re.search(r"Y[ACDEFGHIKLMNPQRSTVWY]{2}L", protein_sequence): # YxxL
        hints.append("YxxL (Tirozin bazlı) motifi: Endositoz sinyali veya hücre içi hedefleme sinyali potansiyeli.")
    if re.search(r"L[ACDEFGHIKLMNPQRSTVWY]L", protein_sequence): # LxL
         hints.append("LxL (Lösin bazlı) motifi: Protein etkileşimleri veya hücre içi hedefleme sinyali potansiyeli (örn: Endoplazmik Retikulum sinyali).")
    if re.search(r"D[ACDEFGHIKLMNPQRSTVWY]{2}L", protein_sequence): # DxxL
         hints.append("DxxL (Aspartat bazlı) motifi: Hücre içi hedefleme sinyali potansiyeli.")
    if re.search(r"LL[ACDEFGHIKLMNPQRSTVWY]{1,4}[ED]", protein_sequence): # Di-leucine
         hints.append("Di-lösin (LL) veya Di-aspartik/glutamik asit (EE/DD) benzeri motifler: Lizozomal/endosomal hedefleme sinyali potansiyeli.")
    # Mitokondriyal Hedefleme Sinyali (Çok Basit - Amfipatik Helix)
    # Yüksek oranda pozitif yüklü (R, K) ve hidrofobik AA'lerin bir arada bulunması
    try: # Bu hesaplamalar hata verebilir
        # isoelectric_point() ve secondary_structure_fraction() analyzer objesi üzerinden çağrılmalı
        if analyzer.isoelectric_point() > 8.0 and analyzer.secondary_structure_fraction()[0] > 0.2 and re.search(r"[KR][KR][ACDEFGHIKLMNPQRSTVWY]{3,9}[AILMFVW]{3,9}", protein_sequence): # Basit amfipatik helix paterni
             hints.append("Yüksek pI, Alfa-Heliks içeriği ve bazik/hidrofobik kümelenme: Mitokondriyal hedefleme sinyali (MTS) potansiyeli.")
    except Exception:
        pass # Eğer pI veya sekonder yapı hesaplanamazsa bu ipucunu atla


    # --- Fizikokimyasal Özelliklere Dayalı İpuçları ---
    # analyzer.instability_index() zaten calculate_protein_stats içinde hesaplanıyor.
    # Değeri direkt oradan alıyoruz.
    try:
        instability_index = analyzer.instability_index()
        if instability_index > 40:
            hints.append(f"Yüksek Kararsızlık İndeksi ({instability_index:.2f}): Kısa ömürlü veya düzenleyici protein olabilir (in vitro).")
        else:
            hints.append(f"Düşük Kararlılık İndeksi ({instability_index:.2f}): Muhtemelen yapısal olarak daha kararlı (in vitro).")
    except Exception as e:
         # hints.append(f"Kararlılık İndeksi hesaplanamadı: {e}") # Hata mesajını ipuçlarına eklemek yerine loglamak daha iyi
         print(f"Kararlılık İndeksi hesaplanamadı: {e}")

    try:
        gravy = analyzer.gravy()
        if gravy < -0.5:
            hints.append(f"Güçlü Hidrofilik (GRAVY={gravy:.3f}): Muhtemelen sitoplazmik, nükleer veya salgılanan çözünür bir protein.")
        elif gravy > 0.5:
            hints.append(f"Güçlü Hidrofobik (GRAVY={gravy:.3f}): Membranla ilişkili veya integral membran proteini olabilir.")
        else:
             hints.append(f"Orta Hidrofobisite (GRAVY={gravy:.3f}): Globüler veya karışık özellikli olabilir.")
    except Exception as e:
         # hints.append(f"Hidrofobisite (GRAVY) hesaplanamadı: {e}")
         print(f"Hidrofobisite (GRAVY) hesaplanamadı: {e}")


    try:
        pI = analyzer.isoelectric_point()
        if pI > 9.0:
            hints.append(f"Yüksek İzoelektrik Nokta (pI={pI:.2f}): Net pozitif yük (nötr pH'da), nükleik asit (DNA/RNA) bağlama potansiyeli.")
        elif pI < 5.0:
            hints.append(f"Düşük İzoelektrik Nokta (pI={pI:.2f}): Net negatif yük (nötr pH'da), metal iyonu bağlama veya asidik proteinlerle etkileşim potansiyeli.")
        else:
             hints.append(f"Nötr İzoelektrik Nokta (pI={pI:.2f}): Nötr veya hafif yüklü (nötr pH'da).")
    except Exception as e:
         # hints.append(f"İzoelektrik Nokta (pI) hesaplanamadı: {e}")
         print(f"İzoelektrik Nokta (pI) hesaplanamadı: {e}")


    # --- Sekonder Yapı Tahminlerine Dayalı İpuçları ---
    try:
        helix, turn, sheet = analyzer.secondary_structure_fraction()
        if helix > 0.4:
            hints.append(f"Yüksek Alfa-Heliks İçeriği (%{helix*100:.1f}): Globüler domainler, transmembran heliksler, DNA bağlama motifleri (örn: helix-turn-helix).")
        if sheet > 0.3:
            hints.append(f"Yüksek Beta-Tabaka İçeriği (%{sheet*100:.1f}): Beta-barrel yapıları (örn: porinler), enzim active bölgeleri, protein etkileşim yüzeyleri.")
        if turn > 0.2: # Turn oranı genelde daha düşüktür
            hints.append(f"Belirgin Dönüş (Turn) İçeriği (%{turn*100:.1f}): Yapısal esneklik, yüzeyde bulunan halkalar (loop).")
    except Exception as e:
        # hints.append(f"İkincil yapı fraksiyonları hesaplanamadı: {e}")
         print(f"İkincil yapı fraksiyonları hesaplanamadı: {e}")


    # --- Genel ---
    # Eğer hiç ipucu bulunamadıysa (sadece başlangıçtaki veri kontrolü veya bir hata mesajı yoksa)
    if len(hints) == 0 or (len(hints) == 1 and ("Analiz için yeterli protein verisi" in hints[0] or "Çok kısa protein dizisi" in hints[0])):
        hints = ["Bu basit analizle belirgin bir fonksiyonel ipucu bulunamadı."]

    # Uyarıyı her zaman sona ekle
    hints.append("UNUTMAYIN: Yukarıdaki ipuçları sadece eğitimsel, spekülatif tahminlerdir ve gerçek fonksiyonel karakterizasyon, deneysel çalışmalar veya ileri düzey biyoinformatik araçlar (BLAST, InterProScan, Pfam vb.) gerektirir.")


    return hints


# --- Ana Analiz Orkestrasyon Fonksiyonu ---
def analyze_sequence_data(sequence_text, options=None):
    """
    Ana analiz fonksiyonu. Sekansı işler, analizleri çalıştırır ve sonuçları toplar.
    """
    if options is None:
        options = {} # Varsayılan opsiyonlar için boş dict

    min_orf_len_aa = options.get('min_orf_length', 50) # JS'den gelen değeri al, yoksa 50
    cpg_algorithm = options.get('cpg_algorithm', 'gardiner') # JS'den gelen algoritma

    results = {
        "original_input": sequence_text[:100] + "..." if len(sequence_text) > 100 else sequence_text, # Girdinin başını göster
        "cleaned_sequence": "",
        "type": "Unknown",
        "length": 0,
        "gc_content": None,
        "distribution": {},
        "start_codons": [], # Pozisyonlar (0-tabanlı)
        "stop_codons": [],  # Pozisyonlar (0-tabanlı)
        "frames": {},       # Çevrilmiş AA dizileri {'+1': 'M...*', ...}
        "orfs": [],         # Bulunan ORF'lar için detaylı dict listesi
        "motifs": [],       # Bulunan motifler için detaylı dict listesi
        "restriction_sites": {}, # {'Enzyme': [pos1, pos2,...]} (0-tabanlı)
        "cpg_islands": [],  # Bulunan CpG adaları için detaylı dict listesi
        "protein_stats": None, # Protein dizisi ise veya en uzun ORF için istatistikler
        "annotations": [],  # Tüm özellikler için birleşik anotasyon listesi (görselleştirme için)
        "analysis_options_used": options, # Hangi opsiyonların kullanıldığını döndür
        "error": None # Genel hata mesajı
    }

    try:
        # 1. Ön İşleme
        cleaned_sequence = preprocess_sequence(sequence_text)
        results["cleaned_sequence"] = cleaned_sequence
        if not cleaned_sequence:
            results["error"] = "Temizleme sonrası sekans boş."
            return results # Boşsa devam etme

        # 2. Sekans Türünü Belirle
        seq_type = detect_sequence_type(cleaned_sequence)
        results["type"] = seq_type

        # 3. Temel İstatistikler
        stats = calculate_stats(cleaned_sequence, seq_type)
        # stats'tan gelen verileri results'a ekle. Hata varsa stats["error"] içinde olacak.
        results["length"] = stats.get("length", 0)
        results["gc_content"] = stats.get("gc_content")
        results["distribution"] = stats.get("distribution", {})
        # stats'tan gelen hatayı genel hata alanına ekle (varsa)
        if stats.get("error"):
            results["error"] = stats.get("error")


        all_annotations = [] # Tüm annotasyonları burada biriktir

        # 4. DNA/RNA Spesifik Analizler
        if seq_type in ['DNA', 'RNA']:
            # 4a. Kodonları Bul ve Çerçeveleri Çevir
            codon_frame_results = find_codons_and_translate_frames(cleaned_sequence, seq_type)
            results["start_codons"] = codon_frame_results.get("start_codons", [])
            results["stop_codons"] = codon_frame_results.get("stop_codons", [])
            results["frames"] = codon_frame_results.get("frames", {})
            all_annotations.extend(codon_frame_results.get("annotations", [])) # Kodon anotasyonlarını ekle

            # 4b. ORF Bul
            # ORF bulma fonksiyonu şimdi doğrudan anotasyon listesi döndürüyor
            orf_annotations = find_orfs(cleaned_sequence, seq_type, min_orf_len_aa)
            # ORF'ları hem kendi listesine hem de genel anotasyon listesine ekle
            results["orfs"] = orf_annotations # Detaylı ORF bilgisi için ayrı liste
            # find_orfs zaten annotation objeleri döndürüyor, bu objeler all_annotations'a eklenecek
            # all_annotations.extend(orf_annotations) # Bu satır aşağıda birleştirme adımında yapılacak


            # 4c. Motif Bul
            # Motif listesi app.py'den options içinde gelebilir veya burada varsayılan kullanılır
            motifs_to_search = options.get('motifs_to_find', None) # options içinden almayı dene
            motif_annotations = find_motifs(cleaned_sequence, seq_type, motifs_to_search)
            results["motifs"] = motif_annotations # Detaylı motif bilgisi
            all_annotations.extend(motif_annotations) # Genel listeye ekle

            # 4d. Restriksiyon Alanları (Sadece DNA)
            if seq_type == 'DNA':
                # Enzim listesi app.py'den options içinde gelebilir
                enzyme_list_to_search = options.get('restriction_enzymes', None) # options içinden almayı dene
                restriction_annotations = find_restriction_sites(cleaned_sequence, seq_type, enzyme_list_to_search)
                # Hata kontrolü ve formatlama
                if restriction_annotations and isinstance(restriction_annotations, list) and restriction_annotations and restriction_annotations[0].get("type") == "error":
                     # Hata durumunda sadece hata mesajını sakla
                     results["restriction_sites"] = {"error": restriction_annotations[0].get("message")}
                     # Hata bilgisini genel hata mesajına ekle
                     results["error"] = (results.get("error", "") + "; " if results.get("error") else "") + restriction_annotations[0].get("message", "Restriksiyon analizi hatası.")

                else:
                    # Başarılı ise formatla: {'Enzyme': [pos1, pos2...]} ve anotasyonları ekle
                    res_sites_dict = {}
                    for anno in restriction_annotations:
                        # Eğer anotasyon bir hata değilse (yukarıdaki if bloğuna girmedik)
                        if anno.get("type") != "error":
                            enzyme = anno.get("enzyme")
                            pos = anno.get("position")
                            if enzyme and pos is not None:
                                 if enzyme not in res_sites_dict:
                                     res_sites_dict[enzyme] = []
                                 # Pozisyonları 0-tabanlı olarak saklıyoruz
                                 res_sites_dict[enzyme].append(pos) # Anotasyondaki pos zaten 0-tabanlı
                            # Anotasyonları genel listeye ekle
                            all_annotations.append(anno) # find_restriction_sites already returns annotation dicts

                    # Pozisyonları sırala (sözlükteki listeler için)
                    for enzyme in res_sites_dict:
                        res_sites_dict[enzyme].sort()
                    results["restriction_sites"] = res_sites_dict


            # 4e. CpG Adaları (Sadece DNA)
            if seq_type == 'DNA':
                # CpG parametrelerini options'tan al
                cpg_options_from_frontend = {
                    'algorithm': cpg_algorithm,
                    'step_size': options.get('cpg_step_size', 1) # Default 1
                }
                # Min length, window size, gc ve oe override'larını ekle
                # None kontrolü yaparak sadece gerçekten override edilenleri geçirelim
                if options.get('cpg_min_length') is not None: cpg_options_from_frontend['min_length'] = options['cpg_min_length']
                if options.get('cpg_window_size') is not None: cpg_options_from_frontend['window_size'] = options['cpg_window_size']
                if options.get('cpg_min_gc') is not None: cpg_options_from_frontend['min_gc_override'] = options['cpg_min_gc']
                if options.get('cpg_min_oe') is not None: cpg_options_from_frontend['min_oe_override'] = options['cpg_min_oe']


                cpg_island_annotations = find_cpg_islands(cleaned_sequence, seq_type, **cpg_options_from_frontend)
                results["cpg_islands"] = cpg_island_annotations # Detaylı ada bilgisi
                all_annotations.extend(cpg_island_annotations) # Genel listeye ekle


        # 5. Protein Spesifik Analizler veya ORF Protein İstatistikleri
        protein_stats_result = None
        if seq_type == 'Protein':
            # Eğer girdi doğrudan protein ise istatistikleri hesapla
            protein_stats_result = calculate_protein_stats(cleaned_sequence) # Bu fonksiyon artık ipuçlarını da içeriyor

        elif seq_type in ['DNA', 'RNA'] and results["orfs"]:
             # En uzun ORF bulunduysa, onun için protein istatistiklerini hesapla
             longest_orf = results["orfs"][0]
             if longest_orf and longest_orf.get("protein_sequence"):
                  protein_stats_result = calculate_protein_stats(longest_orf["protein_sequence"]) # Bu fonksiyon artık ipuçlarını da içeriyor
                  # Protein analizi sırasında hata oluştuysa genel hata alanına ekle
                  if protein_stats_result and protein_stats_result.get("error"):
                       results["error"] = (results.get("error", "") + "; " if results.get("error") else "") + protein_stats_result.get("error")

                  # Protein stats'ı en uzun ORF objesine de ekleyelim (Frontend kolaylığı için)
                  # Dikkat: Bu, yanıtın boyutunu artırabilir.
                  results["orfs"][0]["protein_stats"] = protein_stats_result


        results["protein_stats"] = protein_stats_result # Ana sonuç objesine ekle


        # 6. Tüm Anotasyonları Birleştir ve Sırala
        # all_annotations zaten motifler, kodonlar, cpg adaları, restriksiyon alanlarını içeriyor
        # Şimdi ORF anotasyonlarını da ekleyelim (find_orfs tarafından üretilenler)
        all_annotations.extend(results["orfs"])


        # Anotasyonları pozisyona göre sırala
        # Bazı anotasyonlarda end_position olmayabilir veya start>end olabilir (ters strand)
        # Sıralama için sadece başlangıç pozisyonunu kullanmak genellikle yeterlidir.
        # Veya (start, end) tuple'ına göre sıralanabilir.
        # Anotasyon objelerinin 'position' ve 'end_position' alanları 0-tabanlı olmalı
        # end_position alanı yoksa veya geçerli değilse position ile aynı değeri kullan
        all_annotations.sort(key=lambda x: (x.get('position', -1), x.get('end_position', x.get('position', -1))))

        results["annotations"] = all_annotations


    except Exception as e:
        # analyze_sequence_data fonksiyonunun try bloğunda yakalanmayan herhangi bir hata
        print(f"Beklenmedik Analiz Hatası: {e}\n{traceback.format_exc()}")
        # Var olan hata mesajına bu beklenmedek hatayı ekle
        results["error"] = (results.get("error", "") + "; " if results.get("error") else "") + f"Genel analiz sırasında beklenmedik hata: {e}. Detaylar için sunucu loglarına bakın."

    # Sonuçları döndür
    return results

# Örnek Kullanım (Flask içinde çağrılacak)
if __name__ == '__main__':
    # Test kodları kaldırıldı, burası sadece Flask'i çalıştıracak
    app.run(debug=True, host='0.0.0.0', port=5000)