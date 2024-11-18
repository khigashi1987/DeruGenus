import os
import openai
import json
import time

class LLM():
    def __init__(self, 
                 api_key='', 
                 completion_model_name='gpt-4o-mini'):
        self.client = openai.OpenAI(api_key=api_key)
        self.completion_model_name = completion_model_name

    def openai_wrapper(self,
                       system_setting_prompt='',
                       user_input=''):
        response = self.client.chat.completions.create(
                model=self.completion_model_name,
                messages=[
                    {"role": "system", "content": system_setting_prompt},
                    {"role": "user", "content": user_input},
                ],
                temperature=0.1,
                response_format={ "type": "json_object" },
                seed=8888,
                n=1,
                stop=None,
        )
        try:
            result_json = response.choices[0].message.content.strip()
            if result_json.startswith('```json'):
                result_json = result_json.replace('```json', '')
                result_json = result_json.replace('```', '')
        except Exception as e:
            print(e)
            result_json = None
        return result_json

    def summarize_papers(self,
                         genus_name='',
                         papers_info=[]):
        system_setting_prompt = '''
You are an expert microbiologist with extensive experience in clinical microbiology, microbial ecology, and molecular biology. You have published numerous review papers in major journals and have a deep understanding of both pathogenic and environmental microorganisms. With your expertise, you can effectively evaluate research papers and synthesize key findings about specific microbial taxa.

I have collected abstracts of the top 20 most cited papers where a specific prokaryotic genus name appears. Please analyze these papers and create a concise but thorough review about this microorganism in JSON format, bringing your microbiological expertise to bear while staying grounded in the provided literature.

As a microbiologist, you should:
1. Evaluate the ecological niche and metabolic capabilities of the organism
2. Assess its significance in human health, environmental processes, or biotechnology
3. Consider any unique physiological or genetic characteristics
4. Identify key molecular mechanisms described in the papers
5. Note any methodological advances in studying this organism
6. Recognize emerging research trends in the field

Format your response as follows:

{
  "overview": "A brief description of what kind of microorganism this is, based on the provided papers",
  "research_significance": "Explain its significance - whether it's used as a model organism, causes human/animal/plant diseases, or plays important roles in environmental processes",
  "key_research_topics": "Summarize the main research topics covered in these papers, with reference IDs",
  "research_summary": "A few paragraphs synthesizing the key findings and importance of this organism, citing reference IDs where appropriate",
  "note_on_papers": "If some papers merely mention the genus name but aren't actually about the organism, note which ones to exclude"
}

Please follow these guidelines:
1. Base your summary primarily on the provided abstracts rather than your general knowledge, except for very well-known bacteria like Escherichia, Salmonella, etc.
2. Cite specific Reference IDs when making statements
3. Some papers may just coincidentally mention the genus name - identify and exclude such papers from your analysis
4. Focus on what makes this organism significant in its field
5. Highlight any specific capabilities, metabolic functions, or ecological roles described
6. Note any important applications in research or industry
7. Include relevant technical details but maintain readability for a general scientific audience

The papers are formatted as:

---
Reference ID: [number]
Title: [title]
Citation Count: [count]
Abstract: [abstract]
---

'''

        user_input = f'''
Target Genus: {genus_name}

Papers:
'''
        for i, paper in enumerate(papers_info):
            user_input += f'''
---
Reference ID: {i+1}
Title: {paper['Title']}
Citation Count: {paper['Citation Count']}
Abstract: {paper['Abstract']}
---
'''
        return self.openai_wrapper(system_setting_prompt=system_setting_prompt,
                                   user_input=user_input)

    def translate_info(self,
                       summary_info=[]):
        system_setting_prompt = '''
You are a scientific translator specialized in microbiology and related fields, with extensive experience in translating academic papers and reviews from English to Japanese. You are highly familiar with microbiology terminology in both languages and understand the importance of maintaining technical accuracy while producing natural Japanese text.

Please translate the following JSON content into Japanese. Follow these specific guidelines:

1. Maintain the exact same JSON structure and key names in English
2. Translate only the text content (values) of these keys:
   - "overview"
   - "research_significance"
   - "key_research_topics"
   - "research_summary"
   - "note_on_papers"

3. Important rules for translation:
   - Keep all Reference IDs in their original format (e.g., "Reference IDs: 1, 3, 6" should remain unchanged)
   - Maintain all numerical values as they appear in the original
   - Keep all scientific names in their original Latin form
   - Keep all technical terms consistent throughout the translation
   - Ensure that Japanese text flows naturally while maintaining academic precision
   - Use appropriate Japanese academic writing style (です・ます調)

4. For technical terms, use commonly accepted Japanese translations.

Please provide the complete JSON with translated content while maintaining the exact same structure and formatting as the original.
'''

        summary_info_str = json.dumps(summary_info, indent=4)

        user_input = summary_info_str

        return self.openai_wrapper(system_setting_prompt=system_setting_prompt,
                                   user_input=user_input)




OPENAI_API_KEY = os.environ.get('OPENAI_API_KEY')
COMPLETION_MODEL_NAME = 'gpt-4o-mini'
llmcaller = LLM(api_key=OPENAI_API_KEY,
                completion_model_name=COMPLETION_MODEL_NAME)

TopN = 1000
for ind in range(TopN):
    print(f'Processing {ind}...')

    summary_file_name = f'./summary/summary_{ind}.json'
    if not os.path.exists(summary_file_name):
        json_file = f'./abstracts/abstracts_{ind}.json'
        papers_info = json.load(open(json_file, 'r'))
        genus_name = papers_info['Genus']
        pubmed_hit_count = papers_info['PubMed Hit Count']
        papers = papers_info['Papers']
        # delete paper with empty abstract from papers list
        papers = [paper for paper in papers if paper['Abstract'] != '']
        papers_summary = llmcaller.summarize_papers(genus_name=genus_name,
                                                    papers_info=papers)
        if papers_summary:
            papers_summary = json.loads(papers_summary)
            papers_summary['Genus'] = genus_name
            papers_summary['PubMed Hit Count'] = pubmed_hit_count

            papers_summary['Papers'] = []
            for i, paper in enumerate(papers):
                paper_info = {}
                paper_info['Reference ID'] = i+1
                paper_info['Title'] = paper['Title']
                paper_info['Citation Count'] = paper['Citation Count']
                paper_info['PMID'] = paper['PMID']
                papers_summary['Papers'].append(paper_info)

            with open(summary_file_name, 'w') as f:
                json.dump(papers_summary, f, indent=4)
        else:
            print(f'Failed to summarize.')
        time.sleep(5)
    else:
        print(f'{summary_file_name} already exists.')

    summary_ja_file_name = f'./summary/summary_ja_{ind}.json'
    if not os.path.exists(summary_ja_file_name):
        summary_info = json.load(open(summary_file_name, 'r'))
        summary_info_papers = summary_info['Papers']

        # delete Papers from summary_info
        del summary_info['Papers']

        summary_ja = llmcaller.translate_info(summary_info=summary_info)
        if summary_ja:
            summary_ja = json.loads(summary_ja)
            summary_ja['Papers'] = summary_info_papers
            with open(summary_ja_file_name, 'w', encoding='utf-8') as f:
                json.dump(summary_ja, f, ensure_ascii=False, indent=4)
        else:
            print(f'Failed to translate.')
        time.sleep(5)
    else:
        print(f'{summary_ja_file_name} already exists.')
    