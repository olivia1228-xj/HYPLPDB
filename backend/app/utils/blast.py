import requests


def test_sequence_comparison():
    """Test single sequence comparison"""
    url = 'http://localhost:5000/api/peptides/blast-compare'

    # Test sequence
    test_sequence = "MKWVTFISLLLLFSSAYS"

    # Request payload
    payload = {
        'sequence': test_sequence
    }

    try:
        # Send POST request
        response = requests.post(url, json=payload)

        # Print response status and headers
        print(f"\nSequence Comparison Test:")
        print(f"Status Code: {response.status_code}")
        print(f"Response Headers: {dict(response.headers)}\n")

        # If successful, print results
        if response.status_code == 200:
            results = response.json()
            print("Results:")
            print(f"Query Type: {results.get('types')}")
            print(f"Total Matches: {results.get('total_matches')}")
            print("\nTop 3 Matches:")

            matches = results.get('results', [])
            for match in matches[:3]:  # Print first 3 matches
                print(f"\nID: {match.get('id')}")
                print(f"Name: {match.get('name')}")
                print(f"Score: {match.get('score')}")
                print(f"Identity: {match.get('identity')}%")
        else:
            print(f"Error: {response.text}")

    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")


def test_fasta_comparison():
    """Test FASTA file comparison"""
    url = 'http://localhost:5000/api/peptides/blast-compare'

    # Create a test FASTA file content
    fasta_content = """>Sequence1
MKWVTFISLLLLFSSAYS
>Sequence2
MLCVTVFAAVL"""

    # Create file-like object
    files = {
        'file': ('test.fasta', fasta_content.encode(), 'application/octet-stream')
    }

    try:
        # Send POST request with file
        response = requests.post(url, files=files)

        # Print response status and headers
        print(f"\nFASTA Comparison Test:")
        print(f"Status Code: {response.status_code}")
        print(f"Response Headers: {dict(response.headers)}\n")

        # If successful, print results
        if response.status_code == 200:
            results = response.json()
            print("Results:")
            print(f"Query Type: {results.get('types')}")
            print(f"Total Sequences: {len(results.get('results', []))}")

            # Print details for each sequence
            for idx, sequence_result in enumerate(results.get('results', [])):
                print(f"\nSequence {idx + 1}:")
                print(f"Query Header: {sequence_result.get('query_header')}")
                print(f"Query Sequence: {sequence_result.get('query_sequence')}")
                print("\nTop 2 Matches:")

                matches = sequence_result.get('matches', [])
                for match in matches[:2]:  # Print first 2 matches for each sequence
                    print(f"\nID: {match.get('id')}")
                    print(f"Name: {match.get('name')}")
                    print(f"Score: {match.get('score')}")
                    print(f"Identity: {match.get('identity')}%")
        else:
            print(f"Error: {response.text}")

    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")


if __name__ == "__main__":
    print("Starting BLAST comparison tests...")
    print("=" * 50)

    # Test single sequence comparison
    test_sequence_comparison()

    print("\n" + "=" * 50)

    # Test FASTA file comparison
    test_fasta_comparison()