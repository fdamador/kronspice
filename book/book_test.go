// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package book

import (
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
	"testing"
)

// exampleRef matches a reference to an example folder in a chapter, as a
// relative link (../examples/name) or a path (examples/name).
var exampleRef = regexp.MustCompile(`examples/([A-Za-z0-9_-]+)`)

// TestEveryExampleHasAChapter checks both directions of the book's rule:
// every folder under examples/ is cited by at least one chapter, and every
// example a chapter cites exists.
func TestEveryExampleHasAChapter(t *testing.T) {
	entries, err := os.ReadDir(filepath.Join("..", "examples"))
	if err != nil {
		t.Fatal(err)
	}
	folders := map[string]bool{}
	for _, e := range entries {
		if e.IsDir() {
			folders[e.Name()] = true
		}
	}

	chapters, err := filepath.Glob("[0-9][0-9]-*.md")
	if err != nil {
		t.Fatal(err)
	}
	if len(chapters) == 0 {
		t.Fatal("no chapter files")
	}
	cited := map[string][]string{}
	for _, ch := range chapters {
		text, err := os.ReadFile(ch)
		if err != nil {
			t.Fatal(err)
		}
		for _, m := range exampleRef.FindAllStringSubmatch(string(text), -1) {
			cited[m[1]] = append(cited[m[1]], ch)
		}
	}

	var missing, unknown []string
	for name := range folders {
		if len(cited[name]) == 0 {
			missing = append(missing, name)
		}
	}
	for name, where := range cited {
		if !folders[name] {
			unknown = append(unknown, name+" (in "+strings.Join(where, ", ")+")")
		}
	}
	sort.Strings(missing)
	sort.Strings(unknown)
	if len(missing) > 0 {
		t.Errorf("examples cited by no chapter: %s", strings.Join(missing, ", "))
	}
	if len(unknown) > 0 {
		t.Errorf("chapters cite examples that do not exist: %s", strings.Join(unknown, ", "))
	}
}

// TestContentsListsEveryChapter checks that README.md links every chapter
// file and that every chapter it links exists.
func TestContentsListsEveryChapter(t *testing.T) {
	readme, err := os.ReadFile("README.md")
	if err != nil {
		t.Fatal(err)
	}
	chapters, _ := filepath.Glob("[0-9][0-9]-*.md")
	for _, ch := range chapters {
		if !strings.Contains(string(readme), "("+ch+")") {
			t.Errorf("README.md does not link %s", ch)
		}
	}
	for _, m := range regexp.MustCompile(`\(([0-9][0-9]-[a-z0-9-]+\.md)\)`).FindAllStringSubmatch(string(readme), -1) {
		if _, err := os.Stat(m[1]); err != nil {
			t.Errorf("README.md links %s, which does not exist", m[1])
		}
	}
}
