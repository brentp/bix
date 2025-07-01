package bix

import (
	"testing"
)

func TestExtractEndValue(t *testing.T) {
	tests := []struct {
		name     string
		info     string
		expected int
		hasError bool
	}{
		{
			name:     "END at beginning of INFO",
			info:     "END=12345",
			expected: 12345,
			hasError: false,
		},
		{
			name:     "END at end of INFO",
			info:     "SVTYPE=DEL;SVLEN=-1000;END=12345",
			expected: 12345,
			hasError: false,
		},
		{
			name:     "END in middle of INFO",
			info:     "SVTYPE=DEL;END=12345;SVLEN=-1000",
			expected: 12345,
			hasError: false,
		},
		{
			name:     "END at beginning with multiple fields",
			info:     "END=54321;SVTYPE=DUP;SVLEN=2000",
			expected: 54321,
			hasError: false,
		},
		{
			name:     "No END field",
			info:     "SVTYPE=DEL;SVLEN=-1000;AC=1",
			expected: 0,
			hasError: true,
		},
		{
			name:     "Invalid END value",
			info:     "SVTYPE=DEL;END=not_a_number;AC=1",
			expected: 0,
			hasError: true,
		},
		{
			name:     "Empty INFO field",
			info:     "",
			expected: 0,
			hasError: true,
		},
		{
			name:     "END with zero value",
			info:     "SVTYPE=DEL;END=0",
			expected: 0,
			hasError: false,
		},
		{
			name:     "END field only",
			info:     "END=999",
			expected: 999,
			hasError: false,
		},
		{
			name:     "Multiple semicolons after END",
			info:     "SVTYPE=DEL;END=8888;AC=2;AF=0.5",
			expected: 8888,
			hasError: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := extractEndValue(tt.info)

			if tt.hasError {
				if err == nil {
					t.Errorf("expected error but got none for input: %s", tt.info)
				}
			} else {
				if err != nil {
					t.Errorf("unexpected error for input %s: %v", tt.info, err)
				}
				if result != tt.expected {
					t.Errorf("expected %d but got %d for input: %s", tt.expected, result, tt.info)
				}
			}
		})
	}
}
