import requests


def test_sequence_comparison():
    """Test single sequence comparison"""
    url = 'http://localhost:5000/api/peptides/blast-compare'

    # 修改为使用 form-data 格式
    payload = {
        'sequence': 'MKWVTFISLLLLFSSAYS'
    }

    try:
        # 使用 data 参数而不是 json 参数
        response = requests.post(url, data=payload)

        print(f"\nSequence Comparison Test:")
        print(f"Status Code: {response.status_code}")
        print(f"Response Headers: {dict(response.headers)}\n")

        if response.status_code == 200:
            results = response.json()
            print("Results:")
            for match in results.get('results', [])[:3]:
                print(f"\nID: {match.get('id')}")
                print(f"Species: {match.get('species')}")
                print(f"Score: {match.get('score')}")
                print(f"Identity: {match.get('identity')}")
                print(f"Sequence: {match.get('sequence')}")
        else:
            print(f"Error: {response.text}")

    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")


def test_fasta_comparison():
    """Test FASTA file comparison"""
    url = 'http://localhost:5000/api/peptides/blast-compare'

    # 创建测试 FASTA 内容
    fasta_content = """>Sequence1
MKWVTFISLLLLFSSAYS
>Sequence2
MLCVTVFAAVL"""

    # 创建文件对象
    files = {
        'file': ('test.fasta', fasta_content, 'text/plain')
    }

    try:
        response = requests.post(url, files=files)

        print(f"\nFASTA Comparison Test:")
        print(f"Status Code: {response.status_code}")
        print(f"Response Headers: {dict(response.headers)}\n")

        if response.status_code == 200:
            results = response.json()
            print("Results:")
            print(f"Total Matches: {len(results.get('results', []))}")

            for match in results.get('results', [])[:3]:
                print(f"\nID: {match.get('id')}")
                print(f"Species: {match.get('species')}")
                print(f"Score: {match.get('score')}")
                print(f"Identity: {match.get('identity')}")
                print(f"Sequence: {match.get('sequence')}")
        else:
            print(f"Error: {response.text}")

    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")


if __name__ == "__main__":
    print("Starting BLAST comparison tests...")
    print("=" * 50)

    test_sequence_comparison()
    print("\n" + "=" * 50)
    test_fasta_comparison()