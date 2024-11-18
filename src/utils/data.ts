import genusData from './genus-data.json'

export type Genus = {
  name: string
  overview: string
  research_significance: string
  key_research_topics: string[]
  research_summary: string
  note_on_papers: string
  Papers: {
    "Reference ID": number
    Title: string
    "Citation Count": number
    PMID: string
  }[]
}

export const allGenusData: Genus[] = genusData

export function getAllGenusNames(): string[] {
  return allGenusData.map(genus => genus.name)
}

export function getGenusData(rank: number): Genus | undefined {
  return allGenusData[rank - 1]
}

export function searchGenus(term: string): Genus | null {
  const lowercaseTerm = term.toLowerCase()
  return allGenusData.find(genus => genus.name.toLowerCase() === lowercaseTerm) || null
}