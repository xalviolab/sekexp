# insight_lab/app.py

from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
# import json # Artık doğrudan json kullanmıyoruz

# Analiz fonksiyonlarını sequence_analyzer modülünden import et
from sequence_analyzer import (
    calculate_protein_stats,
    detect_sequence_type,
    find_cpg_islands,
    find_restriction_sites,
    preprocess_sequence,
    find_codons_and_translate_frames,
    find_motifs,
    calculate_stats,
    analyze_sequence_data, # Yeni: Ana orkestrasyon fonksiyonunu import et
    # generate_educational_annotation # Bu fonksiyon artık analyze_sequence_data içinde kullanılıyor
    # find_orfs # Bu fonksiyon artık analyze_sequence_data içinde çağrılıyor
)
# Biopython sınıflarını burada import etmeye gerek yok, sequence_analyzer içinde kalsınlar
# from Bio.Seq import Seq # Gerek yok
# from Bio.SeqUtils import ProtParam # Gerek yok


# --- Flask Uygulama Kurulumu ---
app = Flask(__name__)
CORS(app)

# --- Rotalar ---
@app.route('/')
def index():
    """Ana HTML sayfasını sunar."""
    return render_template('index.html')

@app.route('/analyze-sequence', methods=['POST'])
def analyze_sequence_route():
    """
    Frontend'den gelen sekans verisini alır, analiz eder ve sonuçları JSON olarak döner.
    sequence_analyzer.analyze_sequence_data fonksiyonunu çağırır.
    """
    if not request.is_json:
        return jsonify({"error": "İstek JSON formatında olmalı"}), 400

    data = request.get_json()
    sequence_text = data.get('sequence', None)
    # Frontend'den gelebilecek analiz seçenekleri
    options = data.get('options', {})

    if not sequence_text:
        return jsonify({"error": "Sekans metni sağlanmadı"}), 400

    try:
        # analyze_sequence_data fonksiyonu artık tüm adımları içeriyor
        # Sadece preprocess kontrolü yapıp ana fonksiyona gönderelim
        cleaned_sequence_preview = preprocess_sequence(sequence_text)
        if not cleaned_sequence_preview:
             return jsonify({"error": "Geçerli bir sekans bulunamadı (boş veya sadece başlık içeriyor olabilir)."}), 400

        # Ana analiz fonksiyonunu çağır
        results = analyze_sequence_data(sequence_text, options)

        # analyze_sequence_data kendi içinde hataları yakalayıp 'error' alanına yazıyor
        # Bu nedenle burada response.ok kontrolüne gerek yok, dönen results objesine bakalım
        # if results.get("error"):
            # Hata zaten results objesine eklenmiş durumda
            # return jsonify(results), 500 # Hata varsa 500 döndürmek iyi bir pratik olabilir
            # Ancak frontend'in hata mesajını gösterebilmesi için 200 ile dönmek daha kolay
            # Frontend'in results.error alanını kontrol etmesi yeterli.


        return jsonify(results), 200

    except Exception as e:
        # analyze_sequence_data içinde yakalanmayan beklenmedik bir hata olursa
        print(f"Sunucu hatası: /analyze-sequence endpoint'i çalıştırılırken beklenmedik bir hata oluştu: {e}")
        import traceback
        traceback.print_exc()  # Detaylı hata izini konsola yazdır
        # Genel bir hata mesajı dön
        return jsonify({"error": f"Analiz sırasında beklenmedik bir sunucu hatası oluştu: {e}"}), 500

# --- Uygulamayı Çalıştırma ---
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)