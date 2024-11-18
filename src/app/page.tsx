'use client'

import { useState, useEffect } from 'react'
import { Input } from "@/components/ui/input"
import { Button } from "@/components/ui/button"
import { ScrollArea } from "@/components/ui/scroll-area"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { getAllGenusNames, getGenusData, searchGenus, Genus } from '@/utils/data'

export default function Home() {
  const [selectedGenus, setSelectedGenus] = useState<Genus | null>(null)
  const [selectedRank, setSelectedRank] = useState(0)
  const [genusNames, setGenusNames] = useState<string[]>([])
  const [searchTerm, setSearchTerm] = useState('')
  const [suggestions, setSuggestions] = useState<string[]>([])
  const [activeRange, setActiveRange] = useState('1-100')

  useEffect(() => {
    setGenusNames(getAllGenusNames())
  }, [])

  useEffect(() => {
    if (searchTerm) {
      const filteredSuggestions = genusNames.filter(name =>
        name.toLowerCase().includes(searchTerm.toLowerCase())
      ).slice(0, 10)
      setSuggestions(filteredSuggestions)
    } else {
      setSuggestions([])
    }
  }, [searchTerm, genusNames])

  const handleSearch = (term: string) => {
    const genus = searchGenus(term)
    if (genus) {
      setSelectedGenus(genus)
      setSelectedRank(genusNames.indexOf(genus.name) + 1)
      setActiveRange(getRangeForRank(genusNames.indexOf(genus.name) + 1))
    }
  }

  const handleGenusSelect = (genus: string, rank: number) => {
    const data = getGenusData(rank)
    if (data) {
      setSelectedGenus(data)
      setSelectedRank(rank)
    }
  }

  const getRangeForRank = (rank: number) => {
    return `${Math.floor((rank - 1) / 100) * 100 + 1}-${Math.floor((rank - 1) / 100 + 1) * 100}`
  }

  return (
    <div className="flex h-screen">
      {/* サイドバー */}
      <div className="w-64 bg-gray-100 p-4 overflow-auto">
        <h1 className="text-xl font-bold mb-4">デルジーナス Top1000</h1>
        <div className="mb-4">
          <Input
            type="text"
            placeholder="Genus名を検索..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            className="w-full"
          />
          {suggestions.length > 0 && (
            <ScrollArea className="mt-1 bg-white border border-gray-300 rounded-md shadow-lg max-h-60">
              {suggestions.map((suggestion, index) => (
                <Button
                  key={index}
                  variant="ghost"
                  className="w-full text-left"
                  onClick={() => {
                    setSearchTerm(suggestion)
                    handleSearch(suggestion)
                    setSuggestions([])
                  }}
                >
                  {suggestion}
                </Button>
              ))}
            </ScrollArea>
          )}
        </div>
        <ScrollArea className="h-full">
          {Array.from({ length: 10 }, (_, i) => (
            <Button
              key={i}
              variant={activeRange === `${i * 100 + 1}-${(i + 1) * 100}` ? "default" : "ghost"}
              className="w-full justify-start mb-2"
              onClick={() => setActiveRange(`${i * 100 + 1}-${(i + 1) * 100}`)}
            >
              {`${i * 100 + 1}位〜${(i + 1) * 100}位`}
            </Button>
          ))}
        </ScrollArea>
      </div>

      {/* メインコンテンツエリア */}
      <div className="flex-1 p-4 overflow-auto">
        <ScrollArea className="h-32 mb-4">
          <div className="flex flex-wrap gap-2">
            {genusNames.slice(
              parseInt(activeRange.split('-')[0]) - 1,
              parseInt(activeRange.split('-')[1])
            ).map((genus, index) => (
              <Button
                key={index}
                variant="outline"
                size="sm"
                onClick={() => handleGenusSelect(genus, parseInt(activeRange.split('-')[0]) + index)}
              >
                第{parseInt(activeRange.split('-')[0]) + index}位 {genus}
              </Button>
            ))}
          </div>
        </ScrollArea>

        {selectedGenus && (
          <Card>
            <CardHeader>
              <CardTitle>第{selectedRank}位 {selectedGenus.name}</CardTitle>
              <CardDescription>PubMedヒット回数: {selectedGenus.Papers.reduce((sum, paper) => sum + paper["Citation Count"], 0)}</CardDescription>
            </CardHeader>
            <CardContent>
              <h3 className="font-bold">概要</h3>
              <p>{selectedGenus.overview}</p>
              <h3 className="font-bold mt-2">研究における重要性</h3>
              <p>{selectedGenus.research_significance}</p>
              <h3 className="font-bold mt-2">主要な研究トピック</h3>
              <ul className="list-disc list-inside">
                {selectedGenus.key_research_topics.map((topic, index) => (
                  <li key={index}>{topic}</li>
                ))}
              </ul>
              <h3 className="font-bold mt-2">研究概要</h3>
              <p>{selectedGenus.research_summary}</p>
              <h3 className="font-bold mt-2">参考文献</h3>
              <p className="text-sm text-gray-500">{selectedGenus.note_on_papers}</p>
              <ul className="list-decimal list-inside mt-2">
                {selectedGenus.Papers.map((paper, index) => (
                  <li key={index}>
                    {paper.Title} (引用回数: {paper["Citation Count"]})
                    <a href={`https://pubmed.ncbi.nlm.nih.gov/${paper.PMID}/`} target="_blank" rel="noopener noreferrer" className="text-blue-500 hover:underline ml-2">
                      PubMed
                    </a>
                  </li>
                ))}
              </ul>
            </CardContent>
          </Card>
        )}
      </div>
    </div>
  )
}